import numpy as np
import pytest
from . import *
from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_cartopy_landmask
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift

def test_landmask_cartopy():
    reader_cartopy = reader_cartopy_landmask.Reader ()

    assert reader_cartopy.__on_land__ (np.array([10]), np.array([60])) == [ True ]
    assert reader_cartopy.__on_land__ (np.array([5]), np.array([60])) == [ False]

@pytest.mark.slow
def test_basemap_setup(benchmark):
    benchmark(
        reader_basemap_landmask.Reader,
            llcrnrlon=-1.5, llcrnrlat=59,
            urcrnrlon=7, urcrnrlat=64
        )

def test_cartopy_setup(benchmark):
    benchmark(reader_cartopy_landmask.Reader)

def test_landmask_setup_prepared_cached_interres(benchmark):
    benchmark(reader_cartopy_landmask.Reader, resolution = 'i')

def test_landmask_many(benchmark):
    reader_cartopy = reader_cartopy_landmask.Reader()
    import numpy as np

    y, x = np.mgrid[50:59:50j, -8:2:100j]
    benchmark(reader_cartopy.__on_land__, x, y)

@pytest.mark.slow
def test_cartopy_plot(tmpdir):
    """ Testing cartopy reader against Basemap (directly) """
    reader_cartopy = reader_cartopy_landmask.Reader()

    import matplotlib.pyplot as plt
    import numpy as np
    import cartopy
    import pyproj

    y = np.linspace (20, 80, 360)
    x = np.linspace (0, 60, 360)

    c = reader_cartopy.__on_land__(x, y)

    from mpl_toolkits.basemap import Basemap
    bm = Basemap(projection = 'cyl', resolution = 'c')

    print (reader_cartopy._crs_.globe)

    crs = cartopy.crs.PlateCarree(globe = cartopy.crs.Globe())
    cproj = pyproj.Proj (crs.proj4_init)
    bproj = pyproj.Proj (bm.proj4string)
    print (cproj)
    print (bproj)

    bcrs = cartopy.crs.Mercator(min_latitude = -90, max_latitude=90)
    (bx, by) = pyproj.transform (cproj, bproj, x, y)
    b = [bm.is_land (bbx, bby) for bbx, bby in zip (bx, by)]

    from opendrift_landmask_data import GSHHS
    from shapely import wkb
    with open (GSHHS['f'], 'rb') as fd:
      land = wkb.load(fd)

    fig = plt.figure ()
    ax = plt.axes (projection = cartopy.crs.PlateCarree())
    ax.add_geometries (land, cartopy.crs.PlateCarree(), facecolor = 'none', edgecolor='black')

    ax.scatter (x[c], y[c], c = 'blue', transform = cartopy.crs.PlateCarree(), label = 'cartopy')
    ax.scatter (x[b], y[b], marker = 'x', c = 'green', transform = cartopy.crs.PlateCarree(), label = 'basemap')
    ax.set_global()
    fig.savefig ('%s/cartopy_plot.png' % tmpdir)
    print ("plot saved in %s/cartopy_plot.png" % tmpdir)

    assert len(c) == len(x)
    assert len(c) == len(b)
    assert np.abs((np.isclose(c, b) == False).sum() / len(c)) <= .028


@pytest.mark.slow
def test_basemap_cartopy_matches(test_data):
    reader_cartopy = reader_cartopy_landmask.Reader()
    reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=4, llcrnrlat=59,
                    urcrnrlon=18, urcrnrlat=68.1,
                    resolution='i', projection='merc')

    reader_nordic = reader_ROMS_native.Reader(test_data +
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

    land = (np.array([15.]), np.array([65.6]))
    ocean = (np.array([5.]), np.array([65.6]))

    # basemap
    ob = OceanDrift(loglevel = 00)
    ob.add_reader ([reader_nordic, reader_basemap])

    en, en_prof, missing = ob.get_environment (['land_binary_mask'],
            reader_nordic.start_time,
            land[0], land[1], np.array([0]), None)

    assert en.land_binary_mask == np.array([True])

    en, en_prof, missing = ob.get_environment (['land_binary_mask'],
            reader_nordic.start_time,
            ocean[0], ocean[1], np.array([0]), None)

    assert en.land_binary_mask == np.array([False])
    assert len(ob.readers) == 2

    # cartopy
    oc = OceanDrift(loglevel = 00)
    oc.add_reader ([reader_nordic, reader_cartopy])
    en, en_prof, missing = oc.get_environment (['land_binary_mask'],
            reader_nordic.start_time,
            land[0], land[1], np.array([0]), None)

    assert en.land_binary_mask == np.array([True])

    en, en_prof, missing = oc.get_environment (['land_binary_mask'],
            reader_nordic.start_time,
            ocean[0], ocean[1], np.array([0]), None)

    assert en.land_binary_mask == np.array([False])
    assert len(oc.readers) == 2 # make sure opendrift doesn't add default basemap

def test_cartopy_array(test_data):
    reader_cartopy = reader_cartopy_landmask.Reader()

    reader_nordic = reader_ROMS_native.Reader(test_data +
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

    lon = np.array([15., 5.])
    lat = np.array([65.6, 65.6])

    # cartopy
    oc = OceanDrift(loglevel = 00)
    oc.add_reader ([reader_nordic, reader_cartopy])
    en, en_prof, missing = oc.get_environment (['land_binary_mask'],
            reader_nordic.start_time,
            lon, lat, np.array([0, 0]), None)

    np.testing.assert_array_equal(en.land_binary_mask, np.array([True, False]))
    assert len(oc.readers) == 2 # make sure opendrift doesn't add default basemap

