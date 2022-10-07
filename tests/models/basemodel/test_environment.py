from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.basemodel import *
from opendrift.models.basemodel.environment import Environment

def test_add_reader(test_data):
    reader_arome = reader_netCDF_CF_generic.Reader(test_data / '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
    env = Environment([])
    env.add_reader(reader_arome)

