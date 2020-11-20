import numpy as np
import pytest
import matplotlib.pyplot as plt
from . import *
from opendrift.readers import reader_netCDF_CF_unstructured


def test_open():
    proj = "+proj=utm +zone=33W, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    r = reader_netCDF_CF_unstructured.Reader("niva/AkvaplanNiva_sample.nc",
                                             proj4=proj)
    print(r)
    # r.plot_mesh()
    # plt.show()


def test_contains():
    proj = "+proj=utm +zone=33W, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    r = reader_netCDF_CF_unstructured.Reader("niva/AkvaplanNiva_sample.nc",
                                             proj4=proj)

    x = np.linspace(5.6e5, 6.0e5, 100)
    y = np.linspace(7.76e6, 7.8e6, 100)
    x, y = np.meshgrid(x, y)
    x = x.ravel()
    y = y.ravel()

    # r.plot_mesh()
    # plt.scatter(x, y, marker = 'x', color = 'k')
    # plt.show()

    assert np.all(r.covers_positions(x, y))

    x = np.linspace(5.6e5, 6.0e5, 100)
    y = np.linspace(7.4e6, 7.6e6, 100)
    x, y = np.meshgrid(x, y)
    x = x.ravel()
    y = y.ravel()

    # r.plot_mesh()
    # plt.scatter(x, y, marker = 'x', color = 'b')
    # plt.show()
    assert not np.all(r.covers_positions(x, y))


def test_get_variables():
    proj = "+proj=utm +zone=33W, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    r = reader_netCDF_CF_unstructured.Reader("niva/AkvaplanNiva_sample.nc",
                                             proj4=proj)

    x = np.array([5.6e5])
    y = np.array([7.76e6])
    z = np.array([0])

    u = r.get_variables(['x_sea_water_velocity'], r.start_time, x, y, z)
