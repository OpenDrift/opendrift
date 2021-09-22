import numpy as np
import pytest
import matplotlib.pyplot as plt
from tests import *
from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.models.oceandrift import OceanDrift

akvaplan = "https://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/netcdf_unstructured_samples/AkvaplanNiva_sample_lonlat_fixed.nc"
akvaplan_local = "niva/AkvaplanNiva_sample.nc4"
proj = "+proj=utm +zone=33 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

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

    # test get_variables against manually retrieved value
    tidx = [3,4]
    idx = [150, 170, 171, 2000]
    zidx = [0, 3, 4, 7]

    times = r.times[tidx]
    xc = r.xc[idx]
    yc = r.yc[idx]
    z = (r.dataset['siglay_center'][zidx, idx] * r.dataset['h_center'][idx]).diagonal()
    assert z.shape == xc.shape

    # u is a face / siglay variable
    uu = r.dataset['u'][tidx[0], zidx, idx].diagonal()

    u = r.get_variables(['x_sea_water_velocity'], times[0], xc, yc, z)['x_sea_water_velocity']
    np.testing.assert_array_equal(uu, u)

    uu = r.dataset['u'][tidx[1], zidx, idx].diagonal()

    u = r.get_variables(['x_sea_water_velocity'], times[1], xc, yc, z)['x_sea_water_velocity']
    np.testing.assert_array_equal(uu, u)

def test_get_variables_many(benchmark):
    r = reader_netCDF_CF_unstructured.Reader(akvaplan, proj4=proj)

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

def test_profile():
    o = OceanDrift()
    r = reader_netCDF_CF_unstructured.Reader(akvaplan, proj4=proj)
    o.add_reader(r)
    o.seed_elements(lon=17, lat=70, z=np.atleast_1d([0, -10, -50]),
                    number=3, time=r.start_time)
    print(o)
    o.run(steps=3)
    #o.plot(linecolor='z')

    # Elements at 0 and 10m depth should not have same trajectory
    # This is ok
    assert o.elements.lon[0] != o.elements.lon[1]
    # Elements at 10 and 50m depth should not have same trajectory
    # This is presently failing
    assert o.elements.lon[1] != o.elements.lon[2]
