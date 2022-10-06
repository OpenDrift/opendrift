from .state import State, WrongState
from .init import Init
from .result import Result
from .simulation import Simulation


class OpenDriftSimulation:
    Init = Init
    Simulation = Simulation
    Result = Result

    state: State

    def __init__(self, loglevel=0):
        self.state = Init()

    def run(self):
        match self.state:
            case Init():
                self.state.run()
            case _:
                raise WrongState()

        # TODO: Construct Simulation from Init
        # self.state = Simulation..
        # self.state._run_()
        #
