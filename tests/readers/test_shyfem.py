import numpy as np
import pytest
from tests import *
from opendrift.readers.unstructured import shyfem
from opendrift.models.oceandrift import OceanDrift

ismar = 'https://iws.ismar.cnr.it/thredds/dodsC/emerge/shyfem_unstructured_adriatic.nc'
try:
    r = shyfem.Reader(ismar)
    shyfem_available = True
except:
    shyfem_available = False

#def test_open():
#    r = shyfem.Reader(ismar)
#    print(r)
#    # r.plot_mesh()
#    # plt.show()

need_shyfem = pytest.mark.skipif(shyfem_available == False,
                                 reason='Shyfem sample file not available/readable from iws.ismar.cnr.it/thredds')

@need_shyfem
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

@need_shyfem
def test_z_polarity():
    r = shyfem.Reader(ismar)
    assert r.zmax > r.zmin
    assert (r.z < 0).all()

@need_shyfem
def test_get_values():
    r = shyfem.Reader(ismar)

    x = np.array([16.175186, 16.175186])
    y = np.array([42.457916, 42.457916])
    z = np.array([0, -500])

    u = r.get_variables(['x_sea_water_velocity'], r.start_time, x,
                  y, z)['x_sea_water_velocity']

    ni = 20000
    zi = np.array([0, 25])

    uu = r.dataset['u_velocity'][0, ni, zi]
    np.testing.assert_array_equal(u, uu)

