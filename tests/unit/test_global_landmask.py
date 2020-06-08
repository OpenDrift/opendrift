import numpy as np
import pytest
from . import *
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift

def test_landmask_generate():
    import os, tempfile

    tmpdir = os.path.join (tempfile.gettempdir(), 'landmask')
    mmapf = os.path.join(tmpdir, 'mask.dat')

    if os.path.exists(mmapf): os.unlink(mmapf)

    import opendrift_landmask_data as old
    l = old.Landmask()

    assert os.path.exists(mmapf)

def test_reader_landmask_generate():
    import os, tempfile

    tmpdir = os.path.join (tempfile.gettempdir(), 'landmask')
    mmapf = os.path.join(tmpdir, 'mask.dat')

    if os.path.exists(mmapf): os.unlink(mmapf)

    reader_global = reader_global_landmask.Reader()

    assert os.path.exists(mmapf)

def test_reader_landmask_generate_extent():
    import os, tempfile

    tmpdir = os.path.join (tempfile.gettempdir(), 'landmask')
    mmapf = os.path.join(tmpdir, 'mask.dat')

    if os.path.exists(mmapf): os.unlink(mmapf)

    reader_global = reader_global_landmask.Reader (extent = [4, 55, 11, 65])

    assert os.path.exists(mmapf)

def test_landmask_global():
    reader_global = reader_global_landmask.Reader (extent = [4, 55, 11, 65])

    assert reader_global.__on_land__ (np.array([10]), np.array([60])) == [ True ]
    assert reader_global.__on_land__ (np.array([5]), np.array([60])) == [ False]

def test_global_array(test_data):
    reader_global = reader_global_landmask.Reader()

    reader_nordic = reader_ROMS_native.Reader(test_data +
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

    lon = np.array([15., 5.])
    lat = np.array([65.6, 65.6])

    # global
    oc = OceanDrift(loglevel = 00)
    oc.add_reader ([reader_nordic, reader_global])
    en, en_prof, missing = oc.get_environment (['land_binary_mask'],
            reader_nordic.start_time,
            lon, lat, np.array([0, 0]), None)

    np.testing.assert_array_equal(en.land_binary_mask, np.array([True, False]))
    assert len(oc.readers) == 2 # make sure opendrift doesn't add default basemap


@pytest.mark.veryslow
def test_plot(tmpdir):
    print("setting up global landmask")
    reader_global = reader_global_landmask.Reader(
                        llcrnrlon=18.64, llcrnrlat=69.537,
                        urcrnrlon=19.37, urcrnrlat=69.81)

    x = np.linspace(18.641, 19.369, 10)
    y = np.linspace(69.538, 69.80, 10)

    xx, yy = np.meshgrid(x,y)
    shp = xx.shape
    xx = xx.ravel()
    yy = yy.ravel()

    print ("points:", len(xx))

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy

    reader = cartopy.feature.GSHHSFeature(scale = 'f')

    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    c = reader_global.__on_land__(xx,yy).reshape(shp)
    # c = reader_basemap.__on_land__(xx,yy).reshape(shp)
    print (c)
    ex = [18.641, 19.369, 69.538, 69.80]
    plt.imshow(c, extent = ex, transform = ccrs.PlateCarree())
    ax.coastlines()
    # ax.set_global()
    # plt.show()
    plt.savefig('%s/cartplot.png' % tmpdir)


def test_global_setup(benchmark):
    benchmark(reader_global_landmask.Reader)

@pytest.mark.slow
def test_performance_global(benchmark):
    print("setting up global landmask")
    reader_global = reader_global_landmask.Reader(
            extent = [18.64, 69.537, 19.37, 69.81])


    x = np.linspace(18.641, 19.369, 100)
    y = np.linspace(69.538, 69.80, 100)

    xx, yy = np.meshgrid(x,y)
    xx = xx.ravel()
    yy = yy.ravel()

    print ("points:", len(xx))

    # warmup
    reader_global.__on_land__(xx, yy)

    print ("masking against cartopy")
    benchmark(reader_global.__on_land__, xx,yy)

