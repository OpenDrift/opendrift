from abc import abstractmethod
from .state import State
from .init import Init

class Simulation(State):
    init: Init

    def __init__(self, init: Init):
        self.init = init.copy()

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

    def to_result(self) -> 'Result':
        pass
