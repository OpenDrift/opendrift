import numpy as np
import pytest
import matplotlib.pyplot as plt
from tests import *
from opendrift.readers import reader_netCDF_CF_unstructured

akvaplan = "https://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/netcdf_unstructured_samples/AkvaplanNiva_sample_lonlat_fixed.nc"
akvaplan_local = "niva/AkvaplanNiva_sample.nc"
proj = "+proj=utm +zone=33W, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


def test_open():
    r = reader_netCDF_CF_unstructured.Reader(akvaplan, proj4=proj)
    print(r)
    # r.plot_mesh()
    # plt.show()


def test_contains(benchmark):
    r = reader_netCDF_CF_unstructured.Reader(akvaplan, proj4=proj)

    x = np.linspace(5.6e5, 6.0e5, 100)
    y = np.linspace(7.76e6, 7.8e6, 100)
    x, y = np.meshgrid(x, y)
    x = x.ravel()
    y = y.ravel()

    # r.plot_mesh()
    # plt.scatter(x, y, marker = 'x', color = 'k')
    # plt.show()

    covers = benchmark(r.covers_positions, x, y)
    assert np.all(covers)


def tets_not_contains(benchmark):
    r = reader_netCDF_CF_unstructured.Reader(akvaplan, proj4=proj)

    x = np.linspace(5.6e5, 6.0e5, 100)
    y = np.linspace(7.4e6, 7.6e6, 100)
    x, y = np.meshgrid(x, y)
    x = x.ravel()
    y = y.ravel()

    # r.plot_mesh()
    # plt.scatter(x, y, marker = 'x', color = 'b')
    # plt.show()
    covers = benchmark(r.covers_positions, x, y)
    assert not np.all(covers)


def test_get_variables(benchmark):
    r = reader_netCDF_CF_unstructured.Reader(akvaplan, proj4=proj)

    x = np.array([5.6e5, 5.6e5])
    y = np.array([7.76e6, 7.76e6])
    z = np.array([0, 30])

    u = benchmark(r.get_variables, ['x_sea_water_velocity'], r.start_time, x,
                  y, z)['x_sea_water_velocity']
    print(u)
    print(u.shape)
    assert u.shape == (len(x),)
    assert len(u) == len(x)

@pytest.mark.skip
def test_get_variables_many(benchmark):
    r = reader_netCDF_CF_unstructured.Reader(akvaplan_local, proj4=proj)

    x = np.linspace(5.6e5, 6.0e5, 1000)
    y = np.linspace(7.76e6, 7.8e6, 1000)
    z = np.array([0, 30]).repeat(500)

    u = benchmark(r.get_variables, ['x_sea_water_velocity'], r.start_time, x, y, z)['x_sea_water_velocity']
    assert len(u) == len(x)


def test_nearest_node(benchmark):
    r = reader_netCDF_CF_unstructured.Reader(akvaplan, proj4=proj)
    x = np.linspace(5.6e5, 6.0e5, 100)
    y = np.linspace(7.4e6, 7.6e6, 100)
    x, y = np.meshgrid(x, y)
    x = x.ravel()
    y = y.ravel()
    print("getting nearest for %d points" % len(x))

    benchmark(r._nearest_node_, x, y)


def test_vector_nearest(benchmark):
    X = np.tile(np.arange(0, 30), (50, 1)).T
    xp = np.linspace(0, 30, 50)

    # first dimension is levels, second is positions
    assert X.shape[1] == xp.shape[0]

    nearest = benchmark(reader_netCDF_CF_unstructured.Reader._vector_nearest_,
                        X, xp)

    assert len(nearest) == len(xp)
    truth = [
        0, 1, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 9, 9, 10, 10, 11, 12, 12, 13,
        13, 14, 15, 15, 16, 17, 17, 18, 18, 19, 20, 20, 21, 21, 22, 23, 23, 24,
        24, 25, 26, 26, 27, 28, 28, 29, 29, 29
    ]

    np.testing.assert_array_equal(truth, nearest)
