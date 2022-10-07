from abc import abstractmethod
import logging
from .state import State
from .init import Init

logger = logging.getLogger(__name__)

class Simulation(State):
    ElementType = None

    def __init__(self, init: Init):
        init = init.copy() # TODO: keep this?

        self.ElementType = init.ElementType

    @abstractmethod
    def update(self):
        """
        Step the model one time step.
        """

        return False

    def _run_(self):
        # for ..:
        #   self.update(..)
        pass

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

