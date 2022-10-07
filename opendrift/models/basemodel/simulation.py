from abc import abstractmethod, abstractproperty
from typing import Any, Dict, OrderedDict
import logging
import xarray as xr
import numpy as np
from datetime import timedelta

from opendrift.elements import LagrangianArray

from .environment import Environment
from .state import State
from .init import Init
from .config import Configurable

logger = logging.getLogger(__name__)


class Simulation(State, Configurable):
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


    ## Output
    history: xr.Dataset

    def __init__(self,
                 init: Init,
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
        self.__init_time_steps__(time_step, steps, time_step_output, duration, end_time)

    def __init_time_steps__(self,
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
        if time_step.days < 0:
            logger.info(
                'Backwards simulation, starting from last seeded element')
            self.start_time = self.elements_scheduled_time.max()
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

    @abstractmethod
    def update(self):
        """
        Step the model one time step.
        """

        return False

    def run(self):
        assert len(self.elements_scheduled) > 0, "No elements seeded"

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
