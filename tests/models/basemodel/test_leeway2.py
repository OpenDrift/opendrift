from . import *

from opendrift.models.leeway2 import Leeway

"""Tests for Leeway module."""
def test_leewayprop():
    """Check that Leeway properties are properly read."""
    object_type = 85  # MED-WASTE-7
    lee = Leeway(loglevel=20)
    object_type = object_type
    assert lee.leewayprop[object_type]['Description'] == '>>Medical waste, syringes, small'
    assert lee.leewayprop[object_type]['DWSLOPE'] == 1.79

def test_classes():
    class State:
        pass

    class Model:
        class Init(State):
            def __init__(self):
                print("ModelInit.init")

        class Simulation(State):
            pass

        class Result(State):
            pass

        state: State

        def __init__(self):
            self.state = self.Init()

    class CustomModel(Model):
        class Init(Model.Init):
            def __init__(self):
                super().__init__()
                print("CustomInit.init")

        def __init__(self):
            super().__init__()


    m = Model()

    cm = CustomModel()

