import logging

from .state import State, WrongState
from .init import Init, CONFIG_LEVEL_ADVANCED, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ESSENTIAL
from .result import Result
from .simulation import Simulation

logger = logging.getLogger(__name__)

class OpenDriftSimulation:
    Init = Init
    Simulation = Simulation
    Result = Result

    state: State

    def __init__(self, loglevel=0):
        self.state = self.Init()

    def run(self):
        match self.state:
            case Init():
                self.state = self.Simulation(self.state)

            case Simulation():
                self.state = self.state.run()

            case Result():
                logger.warning("Simulation already done.")

    def __getattr__(self, attr):
        """
        Forward all other method calls and attributes to current state.
        """
        try:
            return getattr(self.state, attr)
        except AttributeError:
            logger.error(f"{attr} does not exist on {self.state}, are you trying to call a method on the wrong state?")
            raise
