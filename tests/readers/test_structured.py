import numpy as np
import xarray as xr
import pandas as pd
import pytest
from . import *
from opendrift.readers import reader_netCDF_CF_generic, reader_ROMS_native
from opendrift.readers.basereader.structured import StructuredReader


svim = 'https://thredds.met.no/thredds/dodsC/nansen-legacy-ocean/SVIM/2020/ocean_avg_20200601.nc4'
try:
    r = xr.open_dataset(svim)
    svim_available = True
except:
    svim_available = False
message = 'SVIM dataset not available from thredds.met.no'


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

@pytest.mark.skipif(svim_available == False, reason=message)
@pytest.mark.slow
def test_lonlat2xy_sequential_big(test_data, benchmark):
    reader = reader_ROMS_native.Reader(svim)
    assert not reader.projected
    reader.__disable_parallel__ = True

    lon = np.arange(-10., 15., .01)
    lat = np.arange(56., 69.8, .01)

    lon, lat = np.meshgrid(lon, lat)
    lon, lat = lon.ravel(), lat.ravel()
    print("coords: %d" % len(lon))

    _, _ = benchmark(reader.lonlat2xy, lon, lat)
    assert reader.__lonlat2xy_parallel__ == False

@pytest.mark.skipif(svim_available == False, reason=message)
@pytest.mark.slow
def test_lonlat2xy_parallel_big(test_data, benchmark):
    reader = reader_ROMS_native.Reader(svim)
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


def test_block_fetched_again_for_deeper_elements():
    # Daily reader with current = depth / 100
    depth = np.array([0., 2, 4, 6, 10, 15, 20, 30, 50])
    u = np.broadcast_to(depth[None, :, None, None] / 100, (2, len(depth), 3, 4)).copy()
    coords = {'time': pd.date_range('2023-01-01', periods=2, freq='1D'),
              'depth': ('depth', depth, {'standard_name': 'depth', 'positive': 'down'}),
              'lat': ('lat', [59., 60, 61], {'standard_name': 'latitude'}),
              'lon': ('lon', [0., 1, 2, 3], {'standard_name': 'longitude'})}
    ds = xr.Dataset({'u': (('time', 'depth', 'lat', 'lon'), u, {'standard_name': 'x_sea_water_velocity'}),
                     'v': (('time', 'depth', 'lat', 'lon'), 0 * u, {'standard_name': 'y_sea_water_velocity'})},
                    coords=coords)
    r = reader_netCDF_CF_generic.Reader(ds)
    variables = ['x_sea_water_velocity', 'y_sea_water_velocity']
    lon, lat = np.array([1., 2.]), np.array([60., 60.])
    time = r.start_time + pd.Timedelta(hours=6)
    for z in ([0, -1], [-30, -1], [0, -1]):
        e, _ = r.get_variables_interpolated(variables, time=time, lon=lon, lat=lat, z=np.array(z, dtype=float))
        np.testing.assert_allclose(e['x_sea_water_velocity'], -np.array(z) / 100, atol=1e-6)
