import numpy as np
import pytest
from . import *
from opendrift.readers import reader_netCDF_CF_generic


def test_covers_positions(test_data):
    reader_arome = reader_netCDF_CF_generic.Reader(
        test_data +
        '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF.nc_20160202_subset')

    ts = reader_arome.get_timeseries_at_position(
        lon=12, lat=68, variables=['x_wind', 'y_wind'])

    assert len(ts['time']) == 49
    x_wind = ts['x_wind']
    assert len(x_wind) == 49

    np.testing.assert_almost_equal(x_wind[0], 2.615, 2)
    np.testing.assert_almost_equal(x_wind[-1], -0.222, 2)


