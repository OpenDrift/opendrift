import numpy as np
import xarray as xr
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
    np.testing.assert_array_almost_equal(x2, np.array([0.082982]))


def test_reader_dataset(test_data):
    ds = xr.open_dataset(test_data + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
    reader_norkyst = reader_netCDF_CF_generic.Reader(ds)
    time = reader_norkyst.start_time
    print(reader_norkyst)


def test_lonlat2xy_sequential(test_data, benchmark):
    reader = reader_ROMS_native.Reader(test_data + '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
    assert not reader.projected
    reader.__disable_parallel__ = True

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

    reader.__disable_parallel__ = True
    xs, ys = reader.lonlat2xy(lon, lat)
    assert reader.__lonlat2xy_parallel__ == False

    reader.__disable_parallel__ = False
    x, y = benchmark(reader.lonlat2xy, lon, lat)
    assert reader.__lonlat2xy_parallel__ == True

    np.testing.assert_equal(x, xs)
    np.testing.assert_equal(y, ys)

@pytest.mark.slow
def test_lonlat2xy_sequential_big(test_data, benchmark):
    reader = reader_ROMS_native.Reader('https://thredds.met.no/thredds/dodsC/nansen-legacy-ocean/SVIM/2020/ocean_avg_20200601.nc4')
    assert not reader.projected
    reader.__disable_parallel__ = True

    lon = np.arange(-10., 15., .01)
    lat = np.arange(56., 69.8, .01)

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = lon.ravel(), lat.ravel()
    print("coords: %d" % len(lon))

    _, _ = benchmark(reader.lonlat2xy, lon, lat)
    assert reader.__lonlat2xy_parallel__ == False

@pytest.mark.slow
def test_lonlat2xy_parallel_big(test_data, benchmark):
    reader = reader_ROMS_native.Reader('https://thredds.met.no/thredds/dodsC/nansen-legacy-ocean/SVIM/2020/ocean_avg_20200601.nc4')
    assert not reader.projected

    lon = np.arange(-10., 15., .01)
    lat = np.arange(56., 69.8, .01)

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = lon.ravel(), lat.ravel()

    reader.__disable_parallel__ = True
    xs, ys = reader.lonlat2xy(lon, lat)
    assert reader.__lonlat2xy_parallel__ == False

    print(xs, ys)
    reader.__disable_parallel__ = False
    x, y = benchmark(reader.lonlat2xy, lon, lat)
    assert reader.__lonlat2xy_parallel__ == True

    np.testing.assert_equal(x, xs)
    np.testing.assert_equal(y, ys)
