from abc import abstractmethod, abstractproperty
from typing import Any, Dict, OrderedDict, List
import types
import traceback
import sys
import logging
import xarray as xr
import numpy as np
from datetime import timedelta, datetime
import pyproj

from opendrift.elements import LagrangianArray
from opendrift.models.physics_methods import PhysicsMethods
from opendrift.timer import Timeable

from .environment import Environment
from .state import State
from .init import Init
from .config import Configurable

logger = logging.getLogger(__name__)


class Simulation(State, Configurable, Timeable, PhysicsMethods):
    ElementType: LagrangianArray

    elements_scheduled: LagrangianArray
    elements_deactivated: LagrangianArray
    elements_scheduled_time: List

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
    stop_on_error = False
    show_continuous_performance = False

    status_categories: List[str]  # Particles are active by default
    minvals: Dict
    maxvals: Dict

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

        iomodule = 'netcdf'
        try:
            io_module = __import__(
                'opendrift.export.io_' + iomodule,
                fromlist=['init', 'write_buffer', 'close', 'import_file'])
        except ImportError:
            logger.info('Could not import iomodule ' + iomodule)
        self.io_init = types.MethodType(io_module.init, self)
        self.io_write_buffer = types.MethodType(io_module.write_buffer, self)
        self.io_close = types.MethodType(io_module.close, self)
        self.io_import_file = types.MethodType(io_module.import_file, self)
        self.io_import_file_xarray = types.MethodType(
            io_module.import_file_xarray, self)

        self.stop_on_error = stop_on_error

        init = init.copy()  # TODO: keep this?
        self._config = init._config

        self.ElementType = init.ElementType

        self.elements_scheduled = init.elements_scheduled
        self.elements_scheduled_time = init.elements_scheduled_time
        self.elements_deactivated = self.ElementType()
        self.elements = self.ElementType()

        if len(self.elements_scheduled) == 0:
            logger.warning("No elements scheduled")

        ## Environment
        self.env = init.env.finalize(self)

        ## Simulation parameters
        self.status_categories = [ "active" ]
        self.minvals = {}  # Dicionaries to store minimum and maximum values of variables
        self.maxvals = {}

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
        self.outfile = outfile
        self.export_variables = export_variables
        self.export_buffer_length = export_buffer_length

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
                                         self.env.required_profiles)

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
                logger.info('========================')
                if self.stop_on_error is True:
                    sys.exit('Stopping on error. ' + self.get_messages())
                if self.steps_calculation <= 1:
                    raise
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
        for var in self.env.required_variables:
            keyword = 'reader_' + var
            if var not in self.env.priority_list:
                if var in self.env.fallback_values:
                    self.add_metadata(keyword, self.env.fallback_values[var])
                else:
                    self.add_metadata(keyword, None)
            else:
                readers = self.env.priority_list[var]
                if readers[0].startswith(
                        'constant_reader') and var in self.env.readers[
                            readers[0]]._parameter_value_map:
                    self.add_metadata(
                        keyword,
                        self.env.readers[readers[0]]._parameter_value_map[var][0])
                else:
                    self.add_metadata(keyword, self.env.priority_list[var])

        if self.outfile is not None:
            logger.debug('Writing and closing output file: %s' % self.outfile)
            # Write buffer to outfile, and close
            if self.steps_output >= self.steps_exported:
                # Write last lines, if needed
                self.io_write_buffer()
            self.io_close()

        # Remove any elements scheduled for deactivation during last step
        self.remove_deactivated_elements()

        if self.export_buffer_length is None:
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
            self.io_import_file(self.outfile)

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

    def store_present_positions(self, IDs=None, lons=None, lats=None):
        """Store present element positions, in case they shall be moved back"""
        if self.get_config('general:coastline_action') == 'previous' or (
                'general:seafloor_action' in self._config
                and self.get_config('general:seafloor_action') == 'previous'):
            if not hasattr(self, 'previous_lon'):
                self.previous_lon = np.ma.masked_all(self.num_elements_total())
                self.previous_lat = np.ma.masked_all(self.num_elements_total())
            if IDs is None:
                IDs = self.elements.ID
                lons = self.elements.lon
                lats = self.elements.lat
                self.newly_seeded_IDs = None
            else:
                # to check if seeded on land
                if len(IDs) > 0:
                    self.newly_seeded_IDs = np.copy(IDs)
                else:
                    self.newly_seeded_IDs = None
            self.previous_lon[IDs - 1] = np.copy(lons)
            self.previous_lat[IDs - 1] = np.copy(lats)

    def increase_age_and_retire(self):
        """Increase age of elements, and retire if older than config setting."""
        # Increase age of elements
        self.elements.age_seconds += self.time_step.total_seconds()

        # Deactivate elements that exceed a certain age
        if self.get_config('drift:max_age_seconds') is not None:
            self.deactivate_elements(self.elements.age_seconds >=
                                     self.get_config('drift:max_age_seconds'),
                                     reason='retired')

        # Deacticate any elements outside validity domain set by user
        if self.validity_domain is not None:
            W, E, S, N = self.validity_domain
            if W is not None:
                self.deactivate_elements(self.elements.lon < W,
                                         reason='outside')
            if E is not None:
                self.deactivate_elements(self.elements.lon > E,
                                         reason='outside')
            if S is not None:
                self.deactivate_elements(self.elements.lat < S,
                                         reason='outside')
            if N is not None:
                self.deactivate_elements(self.elements.lat > N,
                                         reason='outside')

    def interact_with_seafloor(self):
        """Seafloor interaction according to configuration setting"""
        if self.num_elements_active() == 0:
            return
        if 'sea_floor_depth_below_sea_level' not in self.env.priority_list:
            return
        sea_floor_depth = self.sea_floor_depth()
        below = np.where(self.elements.z < -sea_floor_depth)[0]
        if len(below) == 0:
            logger.debug('No elements hit seafloor.')
            return

        i = self.get_config('general:seafloor_action')
        if i == 'lift_to_seafloor':
            logger.debug('Lifting %s elements to seafloor.' % len(below))
            self.elements.z[below] = -sea_floor_depth[below]
        elif i == 'deactivate':
            self.deactivate_elements(self.elements.z < -sea_floor_depth,
                                     reason='seafloor')
            self.elements.z[below] = -sea_floor_depth[below]
        elif i == 'previous':  # Go back to previous position (in water)
            logger.warning('%s elements hit seafloor, '
                           'moving back ' % len(below))
            below_ID = self.elements.ID[below]
            self.elements.lon[below] = \
                np.copy(self.previous_lon[below_ID - 1])
            self.elements.lat[below] = \
                np.copy(self.previous_lat[below_ID - 1])

    def closest_ocean_points(self, lon, lat):
        """Return the closest ocean points for given lon, lat"""

        deltalon = 0.01  # grid
        deltalat = 0.01
        numbuffer = 10
        lonmin = lon.min() - deltalon * numbuffer
        lonmax = lon.max() + deltalon * numbuffer
        latmin = lat.min() - deltalat * numbuffer
        latmax = lat.max() + deltalat * numbuffer
        if not 'land_binary_mask' in self.env.priority_list:
            logger.info('No land reader added, '
                        'making a temporary landmask reader')
            from opendrift.readers import reader_from_url, reader_global_landmask
            land_reader = reader_global_landmask.Reader()
        else:
            logger.info('Using existing reader for land_binary_mask')
            land_reader_name = self.env.priority_list['land_binary_mask'][0]
            land_reader = self.env.readers[land_reader_name]
            o = self

        land = land_reader.__on_land__(lon, lat)

        if land.max() == 0:
            logger.info('All points are in ocean')
            return lon, lat
        logger.info('Moving %i out of %i points from land to water' %
                    (np.sum(land != 0), len(lon)))
        landlons = lon[land != 0]
        landlats = lat[land != 0]
        longrid = np.arange(lonmin, lonmax, deltalon)
        latgrid = np.arange(latmin, latmax, deltalat)
        longrid, latgrid = np.meshgrid(longrid, latgrid)
        longrid = longrid.ravel()
        latgrid = latgrid.ravel()
        # Remove grid-points not covered by this reader
        latgrid_covered = land_reader.covers_positions(longrid, latgrid)[0]
        longrid = longrid[latgrid_covered]
        latgrid = latgrid[latgrid_covered]
        landgrid = o.get_environment(['land_binary_mask'],
                                     lon=longrid,
                                     lat=latgrid,
                                     z=0 * longrid,
                                     time=land_reader.start_time,
                                     profiles=None)[0]['land_binary_mask']
        if landgrid.min() == 1 or np.isnan(landgrid.min()):
            logger.warning('No ocean pixels nearby, cannot move elements.')
            return lon, lat

        oceangridlons = longrid[landgrid == 0]
        oceangridlats = latgrid[landgrid == 0]
        from scipy import spatial
        tree = scipy.spatial.cKDTree(
            np.dstack([oceangridlons, oceangridlats])[0])
        landpoints = np.dstack([landlons, landlats])
        dist, indices = tree.query(landpoints)
        indices = indices.ravel()
        lon[land != 0] = oceangridlons[indices]
        lat[land != 0] = oceangridlats[indices]

        return lon, lat

    def store_previous_variables(self):
        """Store some environment variables, for access at next time step"""

        if not hasattr(self, 'store_previous'):
            return
        if not hasattr(self, 'variables_previous'):
            # Create ndarray to store previous variables
            dtype = [(var, np.float32) for var in self.store_previous]
            self.variables_previous = np.array(np.full(
                self.num_elements_total(), np.nan),
                                               dtype=dtype)

        # Copying variables_previous to environment_previous
        self.environment_previous = self.variables_previous[self.elements.ID -
                                                            1]

        # Use new values for new elements which have no previous value
        for var in self.store_previous:
            undefined = np.isnan(self.environment_previous[var])
            self.environment_previous[var][undefined] = getattr(
                self.environment, var)[undefined]

        self.environment_previous = self.environment_previous.view(np.recarray)

        for var in self.store_previous:
            self.variables_previous[var][self.elements.ID - 1] = getattr(
                self.environment, var)

    def report_missing_variables(self):
        """Issue warning if some environment variables missing."""

        missing_variables = []
        for var in self.env.required_variables:
            if np.isnan(getattr(self.environment, var).min()):
                missing_variables.append(var)

        if len(missing_variables) > 0:
            logger.warning('Missing variables: ' + str(missing_variables))
            self.store_message('Missing variables: ' + str(missing_variables))

    def interact_with_coastline(self, final=False):
        """Coastline interaction according to configuration setting"""
        if self.num_elements_active() == 0:
            return
        i = self.get_config('general:coastline_action')
        if not hasattr(self, 'environment') or not hasattr(
                self.environment, 'land_binary_mask'):
            return
        if i == 'none':  # Do nothing
            return
        if final is True:  # Get land_binary_mask for final location
            en, en_prof, missing = \
                self.env.get_environment(['land_binary_mask'],
                                     self.time,
                                     self.elements.lon,
                                     self.elements.lat,
                                     self.elements.z,
                                     None)
            self.environment.land_binary_mask = en.land_binary_mask

        if i == 'stranding':  # Deactivate elements on land
            self.deactivate_elements(self.environment.land_binary_mask == 1,
                                     reason='stranded')
        elif i == 'previous':  # Go back to previous position (in water)
            if self.newly_seeded_IDs is not None:
                self.deactivate_elements(
                    (self.environment.land_binary_mask == 1) &
                    (self.elements.age_seconds
                     == self.time_step.total_seconds()),
                    reason='seeded_on_land')
            on_land = np.where(self.environment.land_binary_mask == 1)[0]
            if len(on_land) == 0:
                logger.debug('No elements hit coastline.')
            else:
                logger.debug('%s elements hit coastline, '
                             'moving back to water' % len(on_land))
                on_land_ID = self.elements.ID[on_land]
                self.elements.lon[on_land] = \
                    np.copy(self.previous_lon[on_land_ID - 1])
                self.elements.lat[on_land] = \
                    np.copy(self.previous_lat[on_land_ID - 1])
                self.environment.land_binary_mask[on_land] = 0

    def deactivate_elements(self, indices, reason='deactivated'):
        """Schedule deactivated particles for deletion (at end of step)"""
        if any(indices) is False:
            return
        if reason not in self.status_categories:
            self.status_categories.append(reason)
            logger.debug('Added status %s' % (reason))
        reason_number = self.status_categories.index(reason)
        #if not hasattr(self.elements.status, "__len__"):
        if len(np.atleast_1d(self.elements.status)) == 1:
            status = self.elements.status.item()
            self.elements.status = np.zeros(self.num_elements_active())
            self.elements.status.fill(status)
        # Deactivate elements, if they have not already been deactivated
        self.elements.status[indices & (self.elements.status ==0)] = \
            reason_number
        self.elements.moving[indices] = 0
        logger.debug('%s elements scheduled for deactivation (%s)' %
                     (np.sum(indices), reason))
        logger.debug(
            '\t(z: %f to %f)' %
            (self.elements.z[indices].min(), self.elements.z[indices].max()))

    def remove_deactivated_elements(self):
        """Moving deactivated elements from self.elements
        to self.elements_deactivated."""

        # All particles scheduled for deletion
        indices = (self.elements.status != 0)
        #try:
        #    len(indices)
        #except:
        if len(indices) == 0 or np.sum(indices) == 0:
            logger.debug('No elements to deactivate')
            return  # No elements scheduled for deactivation
        # Basic, but some more housekeeping will be required later
        self.elements.move_elements(self.elements_deactivated, indices)
        logger.debug('Removed %i elements.' % (np.sum(indices)))
        if hasattr(self, 'environment'):
            self.environment = self.environment[~indices]
            logger.debug('Removed %i values from environment.' %
                         (np.sum(indices)))
        if hasattr(self, 'environment_profiles') and \
                self.environment_profiles is not None:
            for varname, profiles in self.environment_profiles.items():
                logger.debug('remove items from profile for ' + varname)
                if varname != 'z':
                    self.environment_profiles[varname] = \
                        profiles[:, ~indices]
            logger.debug('Removed %i values from environment_profiles.' %
                         (np.sum(indices)))
            #if self.num_elements_active() == 0:
            #    raise ValueError('No more active elements.')  # End simulation

    def state_to_buffer(self):
        """Append present state (elements and environment) to recarray."""

        steps_calculation_float = \
            (self.steps_calculation * self.time_step.total_seconds() /
             self.time_step_output.total_seconds()) + 1
        if self.time_step <= timedelta(seconds=1):
            self.steps_output = int(np.round(steps_calculation_float))
        else:
            self.steps_output = int(np.floor(steps_calculation_float))

        ID_ind = self.elements.ID - 1
        time_ind = self.steps_output - 1 - self.steps_exported
        if self.steps_calculation == self.expected_steps_calculation:
            final_time_step = True
        else:
            final_time_step = False

        if steps_calculation_float.is_integer() or self.time_step < timedelta(
                seconds=1) or final_time_step is True:
            element_ind = range(len(ID_ind))  # We write all elements
        else:
            deactivated = np.where(self.elements.status != 0)[0]
            if len(deactivated) == 0:
                return  # No deactivated elements this sub-timestep
            # We write history for deactivated elements only:
            logger.debug('Writing history for %s deactivated elements' %
                         len(deactivated))
            ID_ind = ID_ind[deactivated]
            element_ind = deactivated
            time_ind = np.minimum(time_ind + 1, self.history.shape[1] - 1)

        # TODO: storing of variables and environment below should be collected in a single loop
        # Store present state in history recarray
        for i, var in enumerate(self.elements.variables):
            if self.export_variables is not None and \
                    var not in self.export_variables:
                continue
            # Temporarily assuming elements numbered
            # from 0 to num_elements_active()
            # Does not hold when importing ID from a saved file, where
            # some elements have been deactivated
            self.history[var][ID_ind, time_ind] = \
                getattr(self.elements, var)[element_ind]
            if len(ID_ind) > 0:
                newmin = np.min(self.history[var][ID_ind, time_ind])
                newmax = np.max(self.history[var][ID_ind, time_ind])
                if var not in self.minvals:
                    self.minvals[var] = newmin
                    self.maxvals[var] = newmax
                else:
                    self.minvals[var] = np.minimum(self.minvals[var], newmin)
                    self.maxvals[var] = np.maximum(self.maxvals[var], newmax)
        # Copy environment data to history array
        for i, var in enumerate(self.environment.dtype.names):
            if self.export_variables is not None and \
                    var not in self.export_variables:
                continue
            self.history[var][ID_ind, time_ind] = \
                getattr(self.environment, var)[element_ind]
            if len(ID_ind) > 0:
                newmin = np.min(self.history[var][ID_ind, time_ind])
                newmax = np.max(self.history[var][ID_ind, time_ind])
                if var not in self.minvals:
                    self.minvals[var] = newmin
                    self.maxvals[var] = newmax
                else:
                    self.minvals[var] = np.minimum(self.minvals[var], newmin)
                    self.maxvals[var] = np.maximum(self.maxvals[var], newmax)

        # Call writer if buffer is full
        if (self.outfile is not None) and \
                ((self.steps_output - self.steps_exported) ==
                    self.export_buffer_length):
            self.io_write_buffer()


    def update_positions(self, x_vel, y_vel):
        """Move particles according to given velocity components.

        This method shall account for projection metrics (a distance
        on a map projection does not necessarily correspond to the same
        distance over true ground (not yet implemented).

        Arguments:
            x_vel and v_vel: floats, velocities in m/s of particle along
                             x- and y-axes of the inherit SRS (proj4).
        """

        geod = pyproj.Geod(ellps='WGS84')

        azimuth = np.degrees(np.arctan2(x_vel, y_vel))  # Direction of motion
        velocity = np.sqrt(x_vel**2 + y_vel**2)  # Velocity in m/s
        velocity = velocity * self.elements.moving  # Do not move frosen elements

        # Calculate new positions
        self.elements.lon, self.elements.lat, back_az = geod.fwd(
            self.elements.lon, self.elements.lat, azimuth,
            velocity * self.time_step.total_seconds())

        # Check that new positions are valid
        if (self.elements.lon.min() < -180) or (
                self.elements.lon.min() > 360
        ) or (self.elements.lat.min() < -90) or (self.elements.lat.max() > 90):
            logger.info('Invalid new coordinates:')
            logger.info(self.elements)
            raise ValueError()

    def horizontal_diffusion(self):
        """Move elements with random walk according to given horizontal diffuivity."""
        D = self.get_config('drift:horizontal_diffusivity')
        if D == 0:
            logger.debug('Horizontal diffusivity is 0, no random walk.')
            return
        dt = np.abs(self.time_step.total_seconds())
        x_vel = self.elements.moving * np.sqrt(2*D/dt) * np.random.normal(
            scale=1, size=self.num_elements_active())
        y_vel = self.elements.moving * np.sqrt(2*D/dt) * np.random.normal(
            scale=1, size=self.num_elements_active())
        speed = np.sqrt(x_vel * x_vel + y_vel * y_vel)
        logger.debug(
            'Moving elements according to horizontal diffusivity of %s, with speeds between %s and %s m/s'
            % (D, speed.min(), speed.max()))
        self.update_positions(x_vel, y_vel)

