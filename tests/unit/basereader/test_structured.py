import numpy as np
import pytest
from . import *
from opendrift.readers import reader_netCDF_CF_generic

def test_set_convolve(test_data):
    reader_norkyst = reader_netCDF_CF_generic.Reader(test_data + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
    time = reader_norkyst.start_time

    x = reader_norkyst.get_variables_interpolated(['x_sea_water_velocity'], lon = 3., lat = 60., time = time)[0]['x_sea_water_velocity']
    print(x)

    # re-create reader, otherwise input would be cached.
    reader2 = reader_netCDF_CF_generic.Reader(test_data + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
    reader2.set_convolution_kernel(20)
    x2 = reader2.get_variables_interpolated(['x_sea_water_velocity'], lon = 3., lat = 60., time = time)[0]['x_sea_water_velocity']

    assert x != x2
    np.testing.assert_allclose(x2, np.array([0.08291302]))
