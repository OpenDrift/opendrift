from abc import abstractmethod, abstractproperty
from typing import Any, Dict, OrderedDict
import logging
import xarray as xr
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

        ## Environment
        self.env = init.env.finalize(self)


    @abstractmethod
    def update(self):
        """
        Step the model one time step.
        """

        return False

    def run(self,
            time_step=None,
            steps=None,
            time_step_output=None,
            duration=None,
            end_time=None,
            outfile=None,
            export_variables=None,
            export_buffer_length=100,
            stop_on_error=False):
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
