import numpy as np
import pytest
from . import *
from opendrift.readers import reader_netCDF_CF_generic, reader_ROMS_native
from opendrift.readers.basereader.structured import StructuredReader

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

def test_lonlat2xy_sequential(test_data, benchmark):
    reader = reader_ROMS_native.Reader(test_data + '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
    assert not reader.projected
    reader.__parallel_fail__ = True

    lon = np.arange(10., 15., .01)
    lat = np.arange(67.8, 69.8, .01)

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = lon.ravel(), lat.ravel()
    print("coords: %d" % len(lon))

    _, _ = benchmark(reader.lonlat2xy, lon, lat)
    assert reader.__lonlat2xy_parallel__ == False

def test_lonlat2xy_parallel(test_data, benchmark):
    reader = reader_ROMS_native.Reader(test_data + '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
    assert not reader.projected

    lon = np.arange(10., 15., .01)
    lat = np.arange(67.8, 69.8, .01)

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = lon.ravel(), lat.ravel()

    reader.__parallel_fail__ = True
    xs, ys = reader.lonlat2xy(lon, lat)
    assert reader.__lonlat2xy_parallel__ == False

    reader.__parallel_fail__ = False
    x, y = benchmark(reader.lonlat2xy, lon, lat)
    assert reader.__lonlat2xy_parallel__ == True
    assert reader.__parallel_fail__ == False

    np.testing.assert_equal(x, xs)
    np.testing.assert_equal(y, ys)

