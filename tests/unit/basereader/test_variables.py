import numpy as np
import pytest
from . import *
from opendrift.readers import reader_netCDF_CF_generic


# def test_covers_positions(test_data):
#     reader_arome = reader_netCDF_CF_generic.Reader(
#         test_data +
#         '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF.nc_20160202_subset')

#     lon = reader_arome.lon[100, 100]
#     lat = reader_arome.lat[100, 100]
#     ts = reader_arome.get_timeseries_at_position(
#         lon=lon, lat=lat, variables=['x_wind', 'y_wind'])

#     print(reader_arome.times)

#     assert len(ts['time']) == 49
#     x_wind = ts['x_wind']
#     assert len(x_wind) == 49
#     print(x_wind)

#     xx = -1.1851288996433589
#     xx = -1.2184428753105008
#     np.testing.assert_almost_equal(x_wind[0], xx, 2)
#     np.testing.assert_almost_equal(x_wind[-1], -3.991, 2)


# def test_get_variables_interpolated_block(test_data):
#     reader_arome = reader_netCDF_CF_generic.Reader(
#         test_data +
#         '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF.nc_20160202_subset')
#     reader_arome2 = reader_netCDF_CF_generic.Reader(
#         test_data +
#         '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF.nc_20160202_subset')

#     lon = reader_arome.lon[100, 100]
#     lat = reader_arome.lat[100, 100]

#     env, env_profiles = reader_arome.get_variables_interpolated(
#         variables = ['x_wind', 'y_wind'],
#         lon=np.atleast_1d(lon),
#         lat=np.atleast_1d(lat),
#         z = np.atleast_1d(0),
#         time=reader_arome.times[4],
#         rotate_to_proj='+proj=latlon')

#     x_wind = env['x_wind']
#     assert len(x_wind) == 1
#     # np.testing.assert_allclose(x_wind, [-4.54153450])

#     env2, env_profiles2 = reader_arome2.get_variables_interpolated(
#         variables = ['x_wind', 'y_wind'],
#         lon=np.atleast_1d(lon),
#         lat=np.atleast_1d(lat),
#         z = np.atleast_1d(0),
#         time=reader_arome.times[4],
#         rotate_to_proj='+proj=latlon',
#         block = True)

#     x_wind2 = env2['x_wind']

#     np.testing.assert_array_almost_equal(x_wind, x_wind2)
