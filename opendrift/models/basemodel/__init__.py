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

    def run(self, *args, **kwargs):
        match self.state:
            case Init():
                self.into_simulation(*args, **kwargs)
                self.run()

            case Simulation():
                self.state.run()
                self.state = self.Result(self.state)

            case Result():
                logger.warning("Simulation already done.")

    def into_simulation(self, *args, **kwargs):
        match self.state:
            case Init():
                self.state = self.Simulation(self.state, *args, **kwargs)
            case _:
                raise WrongState("Only Init can be converted to Simulation. Try converting to Init first.")

    def __getattr__(self, attr):
        """
        Forward all other method calls and attributes to current state.
        """
        try:
            return getattr(self.state, attr)
        except AttributeError:
            logger.error(f"{attr} does not exist on {self.state}, are you trying to call a method on the wrong state?")
            raise
