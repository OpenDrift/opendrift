from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.basemodel import *

def test_copy_init():
    i = init.Init()
    print(i.copy())

def test_init_hashable():
    i = init.Init()
    print(hash(i))

def test_init_config():
    i = init.Init()
    assert i.get_config('general:time_step_minutes') == 60

def test_construct_simulation():
    i = init.Init()
    s = simulation.Simulation(i, duration = timedelta(days=1))

def test_add_reader(test_data):
    reader_arome = reader_netCDF_CF_generic.Reader(test_data / '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

    i = init.Init()
    i.env.add_reader(reader_arome)

    assert 'arome' in list(i.env.readers.keys())[0]

    s = simulation.Simulation(i, duration = timedelta(days=1))
    assert 'arome' in list(s.env.readers.keys())[0]

    rak = list(s.env.readers.keys())[0]
    ra = s.env.readers[rak]

    assert rak in s.env.readers.keys()
    s.env.discard_reader(ra, 'test')
    assert rak not in s.env.readers.keys()

    # Make sure we didn't mutate the original Environment.
    assert rak in i.env.readers.keys()

