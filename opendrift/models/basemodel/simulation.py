from abc import abstractmethod, abstractproperty
from typing import Any, Dict, OrderedDict
import traceback
import sys
import logging
import xarray as xr
import numpy as np
from datetime import timedelta, datetime

from opendrift.elements import LagrangianArray
from opendrift.timer import Timeable

from .environment import Environment
from .state import State
from .init import Init
from .config import Configurable

logger = logging.getLogger(__name__)


class Simulation(State, Configurable, Timeable):
    ElementType: LagrangianArray

    elements_scheduled: LagrangianArray
    elements_deactivated: LagrangianArray
    """ Active elements being simulated """
    elements: LagrangianArray

    ## Environment
    env: Environment

    ## Simulation parameters
    start_time = None
    time_step = None
    time_step_output = None
    expected_steps_output = None
    expected_steps_calculation = None
    simulation_extent = None

    ## Output
    steps_calculation = 0
    history: xr.Dataset

    def __init__(self,
                 init: Init,
                 start_time=None,
                 time_step=None,
                 steps=None,
                 time_step_output=None,
                 duration=None,
                 end_time=None,
                 outfile=None,
                 export_variables=None,
                 export_buffer_length=100,
                 stop_on_error=False):

        init = init.copy()  # TODO: keep this?
        self._config = init._config

        self.ElementType = init.ElementType

        self.elements_scheduled = init.elements_scheduled
        self.elements_deactivated = self.ElementType()
        self.elements = self.ElementType()

        if len(self.elements_scheduled) == 0:
            logger.warning("No elements scheduled")

        ## Environment
        self.env = init.env.finalize(self)

        ## Simulation parameters
        self.__init_time_steps__(start_time, time_step, steps,
                                 time_step_output, duration, end_time)

        self.expected_steps_output = duration.total_seconds() / \
            self.time_step_output.total_seconds() + 1  # Includes start and end
        self.expected_steps_calculation = duration.total_seconds() / \
            self.time_step.total_seconds()
        self.expected_steps_output = int(self.expected_steps_output)
        self.expected_steps_calculation = int(self.expected_steps_calculation)
        self.expected_end_time = self.start_time + self.expected_steps_calculation * self.time_step

        self.__prepare_simulation_extent__()

        self.__allocate_history_array__(outfile, export_variables, export_buffer_length)

        # Move point seeded on land to ocean
        if self.get_config('seed:ocean_only') is True and \
            ('land_binary_mask' in self.env.required_variables):
            #('land_binary_mask' not in self.fallback_values) and \
            self.timer_start('preparing main loop:moving elements to ocean')
            self.elements_scheduled.lon, self.elements_scheduled.lat = \
                self.closest_ocean_points(self.elements_scheduled.lon,
                                          self.elements_scheduled.lat)
            self.timer_end('preparing main loop:moving elements to ocean')

        #############################
        # Check validity domain
        #############################
        validity_domain = [
            self.get_config('drift:deactivate_west_of'),
            self.get_config('drift:deactivate_east_of'),
            self.get_config('drift:deactivate_south_of'),
            self.get_config('drift:deactivate_north_of')
        ]
        if validity_domain == [None, None, None, None]:
            self.validity_domain = None
        else:
            self.validity_domain = validity_domain

        #############################
        # Model specific preparation
        #############################
        self.prepare_run()

    def __init_time_steps__(self,
                            start_time=None,
                            time_step=None,
                            steps=None,
                            time_step_output=None,
                            duration=None,
                            end_time=None):

        ########################
        # Simulation time step
        ########################
        if time_step is None:
            time_step = timedelta(
                minutes=self.get_config('general:time_step_minutes'))
        if type(time_step) is not timedelta:
            # Time step may be given in seconds, as alternative to timedelta
            time_step = timedelta(seconds=time_step)
        self.time_step = time_step
        if time_step_output is None:
            time_step_output = self.get_config(
                'general:time_step_output_minutes')
            if time_step_output is None:
                self.time_step_output = self.time_step
            else:
                self.time_step_output = timedelta(minutes=time_step_output)
        else:
            if type(time_step_output) is timedelta:
                self.time_step_output = time_step_output
            else:
                self.time_step_output = timedelta(seconds=time_step_output)
            if self.time_step_output.days >= 0 and self.time_step.days < 0:
                self.time_step_output = -self.time_step_output

        time_step_ratio = self.time_step_output.total_seconds() / \
            self.time_step.total_seconds()
        if time_step_ratio < 1:
            raise ValueError('Output time step must be equal or larger '
                             'than calculation time step.')
        if not time_step_ratio.is_integer():
            raise ValueError('Ratio of calculation and output time steps '
                             'must be an integer - given ratio is %s' %
                             time_step_ratio)
        ########################
        # Simulation duration
        ########################
        self.start_time = start_time

        if time_step.days < 0:
            logger.info(
                'Backwards simulation, starting from last seeded element')
            start_time = self.elements_scheduled_time.max()

        if self.start_time != start_time:
            logger.error(
                "Simulation start time does not coincide with first seeded element."
            )

        if (duration is not None and end_time is not None) or \
            (duration is not None and steps is not None) or \
                (steps is not None and end_time is not None):
            raise ValueError('Only one of "steps", "duration" and "end_time" '
                             'may be provided simultaneously')
        if duration is None and end_time is None:
            if steps is not None:
                duration = steps * self.time_step
            else:
                for reader in self.env.readers.values():
                    if reader.end_time is not None:
                        if end_time is None:
                            end_time = reader.end_time
                        else:
                            end_time = min(end_time, reader.end_time)
                    logger.info('Duration, steps or end time not specified, '
                                'running until end of first reader: %s' %
                                (end_time))
        if duration is None:
            duration = end_time - self.start_time

        if time_step.days < 0 and duration.days >= 0:
            # Duration shall also be negative for backwards run
            duration = -duration

        if np.sign(duration.total_seconds()) * np.sign(
                time_step.total_seconds()) < 0:
            raise ValueError(
                "Time step must be negative if duration is negative.")

    def __prepare_simulation_extent__(self):
        ##############################################################
        # Prepare readers for the requested simulation domain/time
        ##############################################################
        max_distance = \
            self.env.max_speed*self.expected_steps_calculation * \
            np.abs(self.time_step.total_seconds())
        deltalat = max_distance / 111000.
        deltalon = deltalat / np.cos(
            np.radians(np.mean(self.elements_scheduled.lat)))
        # TODO: extent should ideally be a general polygon, not only lon/lat-min/max
        # TODO: Should also take into account eventual lifetime of elements
        if len(self.elements_scheduled) > 0:
            simulation_extent = [
                np.maximum(-360,
                           self.elements_scheduled.lon.min() - deltalon),
                np.maximum(-89,
                           self.elements_scheduled.lat.min() - deltalat),
                np.minimum(360,
                           self.elements_scheduled.lon.max() + deltalon),
                np.minimum(89,
                           self.elements_scheduled.lat.max() + deltalat)
            ]
        else:
            simulation_extent = [-360, -89, 360, 89]

        if simulation_extent[2] == 360 and simulation_extent[0] < 0:
            simulation_extent[0] = 0
        logger.debug(
            'Preparing readers for simulation coverage (%s) and time (%s to %s)'
            % (simulation_extent, self.start_time, self.expected_end_time))
        for reader in self.env.readers.values():
            logger.debug('\tPreparing %s' % reader.name)
            reader.prepare(extent=simulation_extent,
                           start_time=self.start_time,
                           end_time=self.expected_end_time,
                           max_speed=self.env.max_speed)
        # Store expected simulation extent, to check if new readers have coverage
        self.simulation_extent = simulation_extent

    def __allocate_history_array__(self,
                                   outfile=None,
                                   export_variables=None,
                                   export_buffer_length=100):
        ####################################################################
        # Preparing history array for storage in memory and eventually file
        ####################################################################
        if export_buffer_length is None:
            self.export_buffer_length = self.expected_steps_output
        else:
            self.export_buffer_length = export_buffer_length

        if self.time_step.days < 0:
            # For backwards simulation, we start at last seeded element
            logger.info('Backwards simulation, starting at '
                        'time of last seeded element')
            self.time = self.elements_scheduled_time.max()
            # Flipping ID array, so that lowest IDs are released first
            self.elements_scheduled.ID = \
                np.flipud(self.elements_scheduled.ID)
        else:
            # Forward simulation, start time has been set when seeding
            self.time = self.start_time

        # Add the output variables which are always required
        if export_variables is not None:
            export_variables = list(
                set(export_variables + ['lon', 'lat', 'ID', 'status']))
        self.export_variables = export_variables
        # Initialise array to hold history (element properties and environment)
        # for export to file.
        history_dtype_fields = [(name,
                                 self.ElementType.variables[name]['dtype'])
                                for name in self.ElementType.variables]
        # Add environment variables
        self.history_metadata = self.ElementType.variables.copy()
        for env_var in self.env.required_variables:
            history_dtype_fields.append((env_var, np.dtype('float32')))
            self.history_metadata[env_var] = {}

        # Remove variables from output array, if only subset is requested
        if self.export_variables is not None:
            history_dtype_fields = [
                f for f in history_dtype_fields
                if f[0] in self.export_variables
            ]
            for m in list(self.history_metadata):
                if m not in self.export_variables:
                    del self.history_metadata[m]

        history_dtype = np.dtype(history_dtype_fields)
        self.history = np.ma.array(np.zeros(
            (len(self.elements_scheduled), self.export_buffer_length)),
                                   dtype=history_dtype)
        self.history.mask = True
        self.steps_exported = 0

        if outfile is not None:
            self.io_init(outfile)
        else:
            self.outfile = None

    @abstractmethod
    def update(self):
        """
        Step the model one time step.
        """

        return False

    @abstractmethod
    def prepare_run(self):
        pass

    def run(self):
        self.add_metadata('simulation_time', datetime.now())
        self.timer_start('main loop')
        for _ in range(self.expected_steps_calculation):
            try:
                # Release elements
                self.release_elements()

                if self.num_elements_active(
                ) == 0 and self.num_elements_scheduled() > 0:
                    self.steps_calculation += 1
                    logger.info(
                        'No active but %s scheduled elements, skipping timestep %s (%s)'
                        % (self.num_elements_scheduled(),
                           self.steps_calculation, self.time))
                    self.state_to_buffer()  # Append status to history array
                    if self.time is not None:
                        self.time = self.time + self.time_step
                    continue

                self.increase_age_and_retire()

                self.interact_with_seafloor()

                if self.show_continuous_performance is True:
                    logger.info(self.performance())
                # Display time to terminal
                logger.debug('===================================' * 2)
                logger.info('%s - step %i of %i - %i active elements '
                            '(%i deactivated)' %
                            (self.time, self.steps_calculation + 1,
                             self.expected_steps_calculation,
                             self.num_elements_active(),
                             self.num_elements_deactivated()))
                logger.debug('%s elements scheduled.' %
                             self.num_elements_scheduled())
                logger.debug('===================================' * 2)

                self.environment, self.environment_profiles, missing = \
                    self.env.get_environment(list(self.env.required_variables),
                                         self.time,
                                         self.elements.lon,
                                         self.elements.lat,
                                         self.elements.z,
                                         self.required_profiles)

                self.store_previous_variables()

                self.calculate_missing_environment_variables()

                if any(missing):
                    self.report_missing_variables()

                self.interact_with_coastline()

                self.interact_with_seafloor()

                self.deactivate_elements(missing, reason='missing_data')

                self.state_to_buffer()  # Append status to history array

                self.remove_deactivated_elements()

                # Propagate one timestep forwards
                self.steps_calculation += 1

                if self.num_elements_active(
                ) == 0 and self.num_elements_scheduled() == 0:
                    raise ValueError(
                        'No more active or scheduled elements, quitting.')

                # Store location, in case elements shall be moved back
                self.store_present_positions()

                #####################################################
                if self.num_elements_active() > 0:
                    logger.debug('Calling %s.update()' % type(self).__name__)
                    self.timer_start('main loop:updating elements')
                    self.update()
                    self.timer_end('main loop:updating elements')
                else:
                    logger.info('No active elements, skipping update() method')
                #####################################################

                self.horizontal_diffusion()

                if self.num_elements_active(
                ) == 0 and self.num_elements_scheduled() == 0:
                    raise ValueError(
                        'No active or scheduled elements, quitting simulation')

                logger.debug('%s active elements (%s deactivated)' %
                             (self.num_elements_active(),
                              self.num_elements_deactivated()))
                # Updating time
                if self.time is not None:
                    self.time = self.time + self.time_step

            except Exception as e:
                message = ('The simulation stopped before requested '
                           'end time was reached.')
                logger.warning(message)
                self.store_message(message)
                logger.info('========================')
                logger.info('End of simulation:')
                logger.info(e)
                logger.info(traceback.format_exc())
                logger.info(self.get_messages())
                if not hasattr(self, 'environment'):
                    sys.exit('Simulation aborted. ' + self.get_messages())
                logger.info('========================')
                if stop_on_error is True:
                    sys.exit('Stopping on error. ' + self.get_messages())
                if self.steps_calculation <= 1:
                    raise ValueError('Simulation stopped within '
                                     'first timestep. ' + self.get_messages())
                break

        self.timer_end('main loop')
        self.timer_start('cleaning up')
        logger.debug('Cleaning up')

        self.interact_with_coastline(final=True)
        self.state_to_buffer()  # Append final status to buffer

        #############################
        # Add some metadata
        #############################
        for var in self.required_variables:
            keyword = 'reader_' + var
            if var not in self.priority_list:
                if var in self.fallback_values:
                    self.add_metadata(keyword, self.fallback_values[var])
                else:
                    self.add_metadata(keyword, None)
            else:
                readers = self.priority_list[var]
                if readers[0].startswith(
                        'constant_reader') and var in self.readers[
                            readers[0]]._parameter_value_map:
                    self.add_metadata(
                        keyword,
                        self.readers[readers[0]]._parameter_value_map[var][0])
                else:
                    self.add_metadata(keyword, self.priority_list[var])

        if outfile is not None:
            logger.debug('Writing and closing output file: %s' % outfile)
            # Write buffer to outfile, and close
            if self.steps_output >= self.steps_exported:
                # Write last lines, if needed
                self.io_write_buffer()
            self.io_close()

        # Remove any elements scheduled for deactivation during last step
        self.remove_deactivated_elements()

        if export_buffer_length is None:
            # Remove columns for unseeded elements in history array
            if self.num_elements_scheduled() > 0:
                logger.info(
                    'Removing %i unseeded elements from history array' %
                    self.num_elements_scheduled())
                mask = np.ones(self.history.shape[0], dtype=bool)
                mask[self.elements_scheduled.ID - 1] = False
                self.history = self.history[mask, :]

            # Remove rows for unreached timsteps in history array
            self.history = self.history[:, range(self.steps_output)]
        else:  # If output has been flushed to file during run, we
            # need to reimport from file to get all data in memory
            del self.environment
            if hasattr(self, 'environment_profiles'):
                del self.environment_profiles
            self.io_import_file(outfile)

        self.timer_end('cleaning up')
        self.timer_end('total time')


    def release_elements(self):
        """Activate elements which are scheduled within following timestep."""

        logger.debug(
            'to be seeded: %s, already seeded %s' %
            (len(self.elements_scheduled), self.num_elements_activated()))
        if len(self.elements_scheduled) == 0:
            return
        if self.time_step.days >= 0:
            indices = (self.elements_scheduled_time >= self.time) & \
                      (self.elements_scheduled_time <
                       self.time + self.time_step)
        else:
            indices = (self.elements_scheduled_time <= self.time) & \
                      (self.elements_scheduled_time >
                       self.time + self.time_step)
        self.store_present_positions(self.elements_scheduled.ID[indices],
                                     self.elements_scheduled.lon[indices],
                                     self.elements_scheduled.lat[indices])
        self.elements_scheduled.move_elements(self.elements, indices)
        self.elements_scheduled_time = self.elements_scheduled_time[~indices]
        logger.debug('Released %i new elements.' % np.sum(indices))

    def num_elements_active(self):
        """The number of active elements."""
        if hasattr(self, 'elements'):
            return len(self.elements)
        else:
            return 0

    def num_elements_deactivated(self):
        """The number of deactivated elements."""
        if hasattr(self, 'elements_deactivated'):
            return len(self.elements_deactivated)
        else:
            return 0

    def num_elements_scheduled(self):
        if hasattr(self, 'elements_scheduled'):
            return len(self.elements_scheduled)
        else:
            return 0

    def num_elements_total(self):
        """The total number of scheduled, active and deactivated elements."""
        return self.num_elements_activated() + self.num_elements_scheduled()

    def num_elements_activated(self):
        """The total number of active and deactivated elements."""
        return self.num_elements_active() + self.num_elements_deactivated()

    def add_metadata(self, key, value):
        """Add item to metadata dictionary, for export as netCDF global attributes"""
        if not hasattr(self, 'metadata_dict'):
            from collections import OrderedDict
            self.metadata_dict = OrderedDict()
        self.metadata_dict[key] = value

    def store_message(self, message):
        """Store important messages to be displayed to user at end."""
        if not hasattr(self, 'messages'):
            self.messages = []
        self.messages.append(message)

    def get_messages(self):
        """Report any messages stored during simulation."""

        if hasattr(self, 'messages'):
            return str(self.messages).strip('[]') + '\n'
        else:
            return ''

