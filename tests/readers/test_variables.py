import numpy as np
import pyproj
from . import *
from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_constant
from opendrift.models.physics_methods import wind_drift_factor_from_trajectory
from opendrift.models.oceandrift import OceanDrift

from opendrift.readers.basereader.variables import ReaderDomain

def test_get_variables_along_trajectory_and_wind_drift_factor_from_trajectory():
    o = OceanDrift(loglevel=50)
    o.add_readers_from_list([o.test_data_folder() +
        '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc',
        o.test_data_folder() +
        '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc'], lazy=False)

    o.env.finalize()
    t = o.env.get_variables_along_trajectory(variables=['x_sea_water_velocity', 'y_sea_water_velocity', 'x_wind', 'y_wind'],
            lons=np.array([3.5, 4, 4.5]), lats=np.array([59.7, 60, 60.3]),
            times=[o.env.readers[list(o.env.readers)[0]].start_time+i*timedelta(hours=3) for i in range(3)])
    np.testing.assert_array_almost_equal(t['x_sea_water_velocity'], [-0.078685, -0.106489, -0.058386])
    np.testing.assert_array_almost_equal(t['x_wind'], [-8.308249, -13.063459, -11.09289])

    wdf, azimuth = wind_drift_factor_from_trajectory(t)
    np.testing.assert_array_almost_equal(wdf, [0.27189012, 0.20492421])
    np.testing.assert_array_almost_equal(azimuth, [73.0112213, 82.39749185])

def test_modulate_longitude_360():
    class R360(ReaderDomain):
        xmin = 0
        xmax = 340
        ymin = -80
        ymax = 80

        def __init__(self):
            self.proj4 = '+proj=lonlat +ellps=WGS84'
            self.crs = pyproj.CRS(self.proj4)
            self.proj = pyproj.Proj(self.proj4)
            super().__init__()

    r = R360()
    lons = np.linspace(0, 300, 100)

    assert (r.modulate_longitude(lons) == lons).all()

    lons = np.array([-180, -90])
    assert (r.modulate_longitude(lons) == np.array([360-180, 360-90])).all()

    lons = np.array([0, 90])
    assert (r.modulate_longitude(lons) == np.array([0, 90])).all()

    lons = np.array([100, 180])
    assert (r.modulate_longitude(lons) == np.array([100, 180])).all()

    lons = np.array([240, 350])
    assert (r.modulate_longitude(lons) == np.array([240, 350])).all()

def test_modulate_longitude_180():
    class R180(ReaderDomain):
        xmin = -150
        xmax = 180
        ymin = -80
        ymax = 80

        def __init__(self):
            self.proj4 = '+proj=lonlat +ellps=WGS84'
            self.crs = pyproj.CRS(self.proj4)
            self.proj = pyproj.Proj(self.proj4)
            super().__init__()

    r = R180()
    lons = np.linspace(-180, 150, 100)

    assert (r.modulate_longitude(lons) == lons).all()

    lons = np.array([-180, -90])
    assert (r.modulate_longitude(lons) == np.array([-180, -90])).all()

    lons = np.array([0, 90])
    assert (r.modulate_longitude(lons) == np.array([0, 90])).all()

    lons = np.array([100, 180])
    assert (r.modulate_longitude(lons) == np.array([100, -180])).all()

    lons = np.array([240])
    assert (r.modulate_longitude(lons) == np.array([-120])).all()


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
    r = reader_constant.Reader({'wind_speed':5, 'wind_to_direction': 225,
                                'land_binary_mask': 0})
    o = OceanDrift(loglevel=50)
    o.set_config('general:use_auto_landmask', False)
    o.add_reader(r)
    o.seed_elements(lon=4, lat=60, time=datetime.now())
    o.run(steps=15)
    np.testing.assert_almost_equal(o.elements.lon, 3.932, 3)
    np.testing.assert_almost_equal(o.elements.lat, 59.966, 3)
    # Wind from SW
    r = reader_constant.Reader({'wind_speed':5, 'wind_to_direction': 45,
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

