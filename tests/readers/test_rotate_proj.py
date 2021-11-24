import pyproj
import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import basereader
from opendrift.models.basemodel import OpenDriftSimulation
from . import *


def test_equiv_proj(test_data):
    """ These two projections are almost equivalent """
    # Only used for projection
    reader_current = reader_netCDF_CF_generic.Reader(
        test_data +
        '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')

    # Simulation SRS
    p = pyproj.Proj(
        '+proj=stere +lat_0=90 +lat_ts=60 +lon_0=70 +x_0=0 +y_0=0 +a=6371000 +rf=298.257223563 +units=m +no_defs'
    )

    lon, lat = 4., 62.
    x, y = reader_current.lonlat2xy(lon, lat)

    # Single value
    u = [0.26724684, 0.2784934]
    ru, rv = reader_current.rotate_vectors(x, y, u[0], u[1],
                                           reader_current.proj, p)

    np.testing.assert_almost_equal(u[0], ru, decimal=3)
    np.testing.assert_almost_equal(u[1], rv, decimal=2)

    # Profile (array of vectors)
    u = np.array([
        0.26724684, 0.26311329, 0.25992924, 0.25909173, 0.25849965, 0.24982239,
        0.25549856, 0.30730924, 0.29996288, 0.23423982, 0.23423982
    ])
    v = np.array([
        0.2784934, 0.26522627, 0.25113547, 0.2446443, 0.23454225, 0.16889542,
        0.18744907, 0.13027644, 0.09087397, 0.19704318, 0.19704318
    ])

    ru, rv = reader_current.rotate_vectors(x, y, u, v, reader_current.proj, p)

    np.testing.assert_array_almost_equal(u, ru, decimal=3)
    np.testing.assert_array_almost_equal(v, rv, decimal=3)


def test_lonlat(test_data):
    # Only used for projection
    reader_current = reader_netCDF_CF_generic.Reader(
        test_data +
        '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')

    # Simulation SRS
    p = OpenDriftSimulation.SRS()

    lon, lat = 4., 62.
    x, y = reader_current.lonlat2xy(lon, lat)

    # Single value
    u = [0.26724684, 0.2784934]
    ru, rv = reader_current.rotate_vectors(x, y, u[0], u[1],
                                           reader_current.proj, p)

    np.testing.assert_almost_equal(-0.14591436614546296, ru, decimal=3)
    np.testing.assert_almost_equal(0.3573351998976778, rv, decimal=2)

    # Profile (array of vectors)
    u = np.array([
        0.26724684, 0.26311329, 0.25992924, 0.25909173, 0.25849965, 0.24982239,
        0.25549856, 0.30730924, 0.29996288, 0.23423982, 0.23423982
    ])
    v = np.array([
        0.2784934, 0.26522627, 0.25113547, 0.2446443, 0.23454225, 0.16889542,
        0.18744907, 0.13027644, 0.09087397, 0.19704318, 0.19704318
    ])

    ru, rv = reader_current.rotate_vectors(x, y, u, v, reader_current.proj, p)


    nu = [
        -0.14591437, -0.13547045, -0.12388817, -0.11829696, -0.10930654,
        -0.05284545, -0.06749336, 0.00579645, 0.0388169, -0.08489617,
        -0.08489617
    ]
    nv = [
        0.3573352, 0.34816854, 0.33953491, 0.33613269, 0.33148788, 0.29689097,
        0.30961478, 0.33373241, 0.31101295, 0.29408664, 0.29408664
    ]
    np.testing.assert_array_almost_equal(nu, ru, decimal=3)
    np.testing.assert_array_almost_equal(nv, rv, decimal=3)

def test_many_rotate_stere(test_data, benchmark):
    # Only used for projection
    reader_current = reader_netCDF_CF_generic.Reader(
        test_data +
        '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')

    # Simulation SRS
    p = OpenDriftSimulation.SRS()

    lon, lat = 4., 62.
    x, y = reader_current.lonlat2xy(lon, lat)

    u = np.linspace(-2, 2, 1000000)
    v = np.linspace(2, -2, 1000000)

    _ru, _rv = benchmark(reader_current.rotate_vectors, x, y, u, v, reader_current.proj, p)

def test_many_rotate_lambert(benchmark):
    # Lambert North
    proj = "+proj=lcc +lat_1=49.50000000000001 +lat_0=49.50000000000001 +lon_0=0 \
            +k_0=0.999877341 +x_0=600000 +y_0=200000 +a=6378249.2 +b=6356515 \
            +units=m +no_defs"

    # Only used for projection
    r = basereader.variables.ReaderDomain()
    r.proj4 = proj
    r.proj = pyproj.Proj(proj)

    # Simulation SRS
    p = OpenDriftSimulation.SRS()

    lon, lat = 4., 62.
    x, y = r.lonlat2xy(lon, lat)

    u = np.linspace(-2, 2, 1000000)
    v = np.linspace(2, -2, 1000000)

    _ru, _rv = benchmark(r.rotate_vectors, x, y, u, v, proj, p)
