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

    assert r.__on_land__(np.array([10]), np.array([60])) == [True]
    assert r.__on_land__(np.array([5]), np.array([60])) == [False]


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
