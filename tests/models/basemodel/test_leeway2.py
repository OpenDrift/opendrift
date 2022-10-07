from opendrift.models.basemodel import *
from opendrift.models.leeway2 import Leeway
from opendrift.readers import reader_netCDF_CF_generic

def test_leeway_prop():
    """Check that Leeway properties are properly read."""
    object_type = 85  # MED-WASTE-7
    lee = Leeway(loglevel=20)
    object_type = object_type

    assert isinstance(lee.state, Init)
    assert lee.state.leewayprop[object_type]['Description'] == '>>Medical waste, syringes, small'
    assert lee.state.leewayprop[object_type]['DWSLOPE'] == 1.79

def test_leeway_init():
    l = Leeway()

def test_leeway_seed_elements(test_data):
    l = Leeway()
    reader_arome = reader_netCDF_CF_generic.Reader(test_data / '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
    l.seed_elements(lon=4.5, lat=59.6, radius=100, number=1000,
                    time=reader_arome.start_time, object_type=26)

def test_leeway_construct_simulation(test_data):
    l = Leeway()
    reader_arome = reader_netCDF_CF_generic.Reader(test_data / '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
    l.seed_elements(lon=4.5, lat=59.6, radius=100, number=1000,
                    time=reader_arome.start_time, object_type=26)

    l.into_simulation()

def test_leeway_run_simulation(test_data):
    l = Leeway()
    reader_arome = reader_netCDF_CF_generic.Reader(test_data / '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
    l.seed_elements(lon=4.5, lat=59.6, radius=100, number=1000,
                    time=reader_arome.start_time, object_type=26)

    l.into_simulation()
    l.run()

