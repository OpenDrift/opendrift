from .init import Init
from .result import Result
from .simulation import Simulation

class OpenDriftSimulation:
    Init = Init
    state = None

    def __init__(self, loglevel=0):
        self = Init()

