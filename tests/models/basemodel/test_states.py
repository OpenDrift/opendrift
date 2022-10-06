from opendrift.models import basemodel

def test_copy_init():
    i = basemodel.init.Init()
    print(i.copy())

def test_init_hashable():
    i = basemodel.init.Init()
    print(hash(i))
