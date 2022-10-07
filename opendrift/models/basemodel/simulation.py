from abc import abstractmethod, abstractproperty
from typing import Any
import logging
import xarray as xr
from opendrift.elements import LagrangianArray
from .state import State
from .init import Init

logger = logging.getLogger(__name__)

class Simulation(State):
    ElementType: LagrangianArray

    elements_scheduled: LagrangianArray
    elements_deactivated: LagrangianArray

    """ Active elements being simulated """
    elements: LagrangianArray

    ## Environment
    fallback_values: dict

    ## Output
    history: xr.Dataset

    def __init__(self, init: Init):
        init = init.copy() # TODO: keep this?

        self.ElementType = init.ElementType

        self.elements_scheduled = init.elements_scheduled
        self.elements_deactivated = self.ElementType()
        self.elements = self.ElementType()

        self.fallback_values = init.get_fallback_values()



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

