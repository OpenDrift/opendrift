import pyproj
import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from . import *

def test_rotate_vector (test_data):
    """ These two projections are almost equivalent """
    # Only used for projection
    reader_current = reader_netCDF_CF_generic.Reader(test_data +
            '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')

    # Simulation SRS
    p = pyproj.Proj ('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=70 +x_0=0 +y_0=0 +a=6371000 +rf=298.257223563 +units=m +no_defs')

    lon, lat = 4., 62.
    x, y = reader_current.lonlat2xy (lon, lat)

    # Single value
    u = [0.26724684, 0.2784934]
    ru, rv = reader_current.rotate_vectors (
            x, y,
            u[0], u[1],
            reader_current.proj,
            p)

    np.testing.assert_almost_equal (u[0], ru, decimal=3)
    np.testing.assert_almost_equal (u[1], rv, decimal=2)

    # Profile (array of vectors)
    u = np.array([0.26724684, 0.26311329, 0.25992924, 0.25909173, 0.25849965, 0.24982239, 0.25549856, 0.30730924, 0.29996288, 0.23423982, 0.23423982])
    v = np.array([0.2784934 , 0.26522627, 0.25113547, 0.2446443 , 0.23454225, 0.16889542, 0.18744907, 0.13027644, 0.09087397, 0.19704318, 0.19704318])

    ru, rv = reader_current.rotate_vectors (
            x, y,
            u, v,
            reader_current.proj,
            p)

    np.testing.assert_array_almost_equal (u, ru, decimal = 3)
    np.testing.assert_array_almost_equal (v, rv, decimal = 3)


