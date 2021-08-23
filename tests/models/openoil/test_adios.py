from opendrift.models.openoil import adios

def test_get_all():
    oils = adios.oils(1)
    print(oils)
