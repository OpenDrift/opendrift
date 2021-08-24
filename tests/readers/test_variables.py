import numpy as np
import pytest
from . import *
from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_constant
from opendrift.models.oceandrift import OceanDrift


def test_covers_positions(test_data):
    reader_arome = reader_netCDF_CF_generic.Reader(
        test_data +
        '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF_20160202_subset.nc')

    ts = reader_arome.get_timeseries_at_position(
        lon=12, lat=68, variables=['x_wind', 'y_wind'])

    assert len(ts['time']) == 49
    x_wind = ts['x_wind']
    assert len(x_wind) == 49

    np.testing.assert_almost_equal(x_wind[0], 2.836, 2)
    np.testing.assert_almost_equal(x_wind[-1], -0.667, 2)

def test_environment_mapping(test_data):

    # Wind from NE
    r = reader_constant.Reader({'wind_speed':5, 'wind_from_direction': 45,
                                'land_binary_mask': 0})
    o = OceanDrift(loglevel=50)
    o.set_config('general:use_auto_landmask', False)
    o.add_reader(r)
    o.seed_elements(lon=4, lat=60, time=datetime.now())
    o.run(steps=15)
    np.testing.assert_almost_equal(o.elements.lon, 3.932, 3)
    np.testing.assert_almost_equal(o.elements.lat, 59.966, 3)
    # Wind from SW
    r = reader_constant.Reader({'wind_speed':5, 'wind_from_direction': 225,
                                'land_binary_mask': 0})
    o = OceanDrift(loglevel=50)
    o.set_config('general:use_auto_landmask', False)
    o.add_reader(r)
    o.seed_elements(lon=4, lat=60, time=datetime.now())
    o.run(steps=15)
    np.testing.assert_almost_equal(o.elements.lon, 4.068, 3)
    np.testing.assert_almost_equal(o.elements.lat, 60.034, 3)

    # land_binary_mask mapped from sea_floor_depth_below_sea_level
    r = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
            '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
    assert 'land_binary_mask' not in r.derived_variables  # Disabled by default
    r.activate_environment_mapping('land_binary_mask_from_ocean_depth')
    assert 'land_binary_mask' in r.derived_variables

