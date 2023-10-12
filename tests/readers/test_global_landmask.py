import numpy as np
import pytest
from . import *
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift

def test_global_setup(benchmark):
    benchmark(reader_global_landmask.Reader)

def test_landmask_global():
    reader_global = reader_global_landmask.Reader()

    assert reader_global.__on_land__(np.array([10]), np.array([60])) == [True]
    assert reader_global.__on_land__(np.array([5]), np.array([60])) == [False]


def test_global_array(test_data):
    reader_global = reader_global_landmask.Reader()

    reader_nordic = reader_ROMS_native.Reader(
        test_data +
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

    lon = np.array([15., 5.])
    lat = np.array([65.6, 65.6])

    # global
    oc = OceanDrift(loglevel=00)
    oc.add_reader([reader_nordic, reader_global])
    oc.env.finalize()
    en, en_prof, missing = oc.env.get_environment(['land_binary_mask'],
                                              reader_nordic.start_time, lon,
                                              lat, np.array([0, 0]), None)

    np.testing.assert_array_equal(en.land_binary_mask, np.array([True, False]))
    assert len(
        oc.env.readers) == 2  # make sure opendrift doesn't add default basemap


@pytest.mark.veryslow
def test_plot(tmpdir):
    print("setting up global landmask")
    reader_global = reader_global_landmask.Reader()

    x = np.linspace(18.641, 19.369, 10)
    y = np.linspace(69.538, 69.80, 10)

    xx, yy = np.meshgrid(x, y)
    shp = xx.shape
    xx = xx.ravel()
    yy = yy.ravel()

    print("points:", len(xx))

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    c = reader_global.__on_land__(xx, yy).reshape(shp)
    # c = reader_basemap.__on_land__(xx,yy).reshape(shp)
    print(c)
    ex = [18.641, 19.369, 69.538, 69.80]
    plt.imshow(c, extent=ex, transform=ccrs.PlateCarree())
    ax.coastlines()
    # ax.set_global()
    # plt.show()
    plt.savefig('%s/cartplot.png' % tmpdir)

@pytest.mark.slow
@pytest.mark.parametrize("fast", [False, True])
@pytest.mark.parametrize("scale", ["auto", "c", "f"])
@pytest.mark.mpl_image_compare(remove_text=True, tolerance=2.35)
def test_plot_auto_scale(test_data, scale, fast, show_plot):
    reader_global = reader_global_landmask.Reader()
    reader_nordic = reader_ROMS_native.Reader(
        test_data +
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

    oc = OceanDrift(loglevel=0)
    oc.add_reader([reader_nordic, reader_global])
    oc.seed_elements(lon=4.8,
                    lat=60.0,
                    number=10,
                    radius=1000,
                    time=reader_nordic.start_time)
    oc.run(steps=2)

    fig = oc.plot(buffer=5., lscale=scale, fast=fast, show=show_plot)[1]
    print(fig)
    return fig

@pytest.mark.slow
def test_performance_global(benchmark):
    print("setting up global landmask")
    reader_global = reader_global_landmask.Reader()

    x = np.linspace(18.641, 19.369, 100)
    y = np.linspace(69.538, 69.80, 100)

    xx, yy = np.meshgrid(x, y)
    xx = xx.ravel()
    yy = yy.ravel()

    print("points:", len(xx))

    # warmup
    reader_global.__on_land__(xx, yy)

    benchmark(reader_global.__on_land__, xx, yy)

def test_dateline():
    mask = reader_global_landmask.Reader()

    x = np.linspace(-180, 180, 100)
    y = np.linspace(-90, 90, 100)

    xx, yy = np.meshgrid(x, y)
    xx, yy = xx.ravel(), yy.ravel()
    mm = mask.__on_land__(xx, yy)

    # Offset
    x2 = np.linspace(180, 540, 100)
    y2 = np.linspace(-90, 90, 100)


    xx, yy = np.meshgrid(x2, y2)
    xx, yy = xx.ravel(), yy.ravel()
    MM = mask.__on_land__(xx, yy)

    np.testing.assert_array_equal(mm, MM)
