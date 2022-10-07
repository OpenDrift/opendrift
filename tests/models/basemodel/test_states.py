from opendrift.models.basemodel import *

def test_copy_init():
    i = init.Init()
    print(i.copy())

def test_init_hashable():
    i = init.Init()
    print(hash(i))

def test_construct_simulation():
    i = init.Init()
    s = simulation.Simulation(i)
