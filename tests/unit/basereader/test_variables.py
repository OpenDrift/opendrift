import numpy as np
import pytest
from . import *
from opendrift.readers import reader_netCDF_CF_generic


def test_covers_positions(test_data):
    reader_arome = reader_netCDF_CF_generic.Reader(
        test_data +
        '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF.nc_20160202_subset')

    lon = reader_arome.lon[100, 100]
    lat = reader_arome.lat[100, 100]
    ts = reader_arome.get_timeseries_at_position(
        lon=lon, lat=lat, variables=['x_wind', 'y_wind'])

    assert len(ts['time']) == 49
    x_wind = ts['x_wind']
    assert len(x_wind) == 49

    np.testing.assert_almost_equal(x_wind[0], -1.185, 2)
    np.testing.assert_almost_equal(x_wind[-1], -3.870, 2)


