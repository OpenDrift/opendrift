import numpy as np
import pytest
from . import *
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_shape
import cartopy.io.shapereader as shpreader
import urllib.error

has_ne_shapes = False
try:
    shpreader.natural_earth(resolution='110m',
                            category='cultural',
                            name='admin_0_countries')
    has_ne_shapes = True
except urllib.error.URLError:
    has_ne_shapes = False

need_ne_shapes = pytest.mark.skipif(has_ne_shapes == False,
                                    reason='Natural Earth shapes are not available.')


@need_ne_shapes
def test_on_land():
    shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')
    r = reader_shape.Reader.from_shpfiles(shpfilename)

    assert r.__on_land__(np.array([10]), np.array([60]), r.variables[0]) == [True]
    assert r.__on_land__(np.array([5]), np.array([60]), r.variables[0]) == [False]


@need_ne_shapes
def test_global_array(test_data):
    shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')
    reader_landmask = reader_shape.Reader.from_shpfiles(shpfilename)

    reader_nordic = reader_ROMS_native.Reader(
        test_data +
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

    lon = np.array([15., 5.])
    lat = np.array([65.6, 65.6])

    # global
    oc = OceanDrift(loglevel=00)
    oc.set_config('general:use_auto_landmask', False)
    oc.add_reader([reader_nordic, reader_landmask])
    oc.env.finalize()
    en, en_prof, missing = oc.env.get_environment(['land_binary_mask'],
                                              reader_nordic.start_time, lon,
                                              lat, np.array([0, 0]), None)

    np.testing.assert_array_equal(en.land_binary_mask, np.array([True, False]))
    print(oc.env.readers)
    assert len(
        oc.env.readers) == 2  # make sure opendrift doesn't add default basemap

def test_get_shape_variable():
    oc = OceanDrift(loglevel=30)
    shpfilename = [
        oc.test_data_folder() + 'shp/test_vars.shp',
        oc.test_data_folder() + 'shp/test_more_vars.shp',
    ]
    reader = reader_shape.Reader.from_shpfiles(shpfilename)
    test = reader.get_variables(['var0', 'var1'],
        x=14.267009, y=68.528457)
    assert test['var0'] == False and test['var1'] == True, \
        "Failed on reading test_vars.shp"
    test = reader.get_variables(['var0', 'var1'],
        x=13.979286, y=68.164405)
    assert test['var0'] == True and test['var1'] == False, \
        "Failed on reading test_more_vars.shp"

def test_global_var_array():
    oc = OceanDrift(loglevel=30)
    shpfilename = [
        oc.test_data_folder() + 'shp/test_vars.shp',
        oc.test_data_folder() + 'shp/test_vars_landmask.shp',
        oc.test_data_folder() + 'shp/test_more_vars.shp',
    ]
    reader_landmask = reader_shape.Reader.from_shpfiles(shpfilename)

    reader_nordic = reader_ROMS_native.Reader(
        oc.test_data_folder() + \
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

    lon = np.array([15., 5.])
    lat = np.array([65.6, 65.6])

    # global
    oc = OceanDrift(loglevel=00)
    oc.set_config('general:use_auto_landmask', False)
    oc.add_reader([reader_nordic, reader_landmask])
    oc.env.finalize()
    en, en_prof, missing = oc.env.get_environment(['land_binary_mask'],
                                              reader_nordic.start_time, lon,
                                              lat, np.array([0, 0]), None)

    np.testing.assert_array_equal(en.land_binary_mask, np.array([True, False]))
    print(oc.env.readers)
    assert len(
        oc.env.readers) == 2  # make sure opendrift doesn't add default basemap

