import numpy as np
import pytest
from tests import *
from opendrift.readers.unstructured import shyfem
from opendrift.models.oceandrift import OceanDrift

ismar = 'https://iws.ismar.cnr.it/thredds/dodsC/emerge/shyfem_unstrct_ADRIA_20210408.nc'


def test_open():
    r = shyfem.Reader(ismar)
    print(r)
    # r.plot_mesh()
    # plt.show()


def test_get_variables(benchmark):
    r = shyfem.Reader(ismar)

    x = np.array([13., 13.])
    y = np.array([40, 40.1])
    z = np.array([0, -100])

    u = benchmark(r.get_variables, ['x_sea_water_velocity'], r.start_time, x,
                  y, z)['x_sea_water_velocity']
    print(u)
    print(u.shape)
    assert u.shape == (len(x),)
    assert len(u) == len(x)

def test_z_polarity():
    r = shyfem.Reader(ismar)
    assert r.zmax > r.zmin
    assert (r.z < 0).all()

def test_get_values():
    r = shyfem.Reader(ismar)

    x = np.array([12.411164, 12.411164])
    y = np.array([45.445023, 45.445023])
<<<<<<< HEAD
    z = np.array([-1.6, -500])
=======
    z = np.array([0, -500])
>>>>>>> gauteh/shyfem

    u = r.get_variables(['x_sea_water_velocity'], r.start_time, x,
                  y, z)['x_sea_water_velocity']

    ni = 20000
<<<<<<< HEAD
    zi = np.array([0, -25])
=======
    zi = np.array([0, 25])
>>>>>>> gauteh/shyfem

    uu = r.dataset['u_velocity'][0, ni, zi]
    np.testing.assert_array_equal(u, uu)

