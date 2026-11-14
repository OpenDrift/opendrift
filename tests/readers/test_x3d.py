from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic

pytest.importorskip('x3d')

x3d_url = ('x3d+https://rgw.met.no/5e4f32366fc74342866b4c67dbd4713d:'
           'norkyst800m-v3/norkyst800')
thredds_url = 'https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be'

try:
    ds_x3d = xr.open_dataset(x3d_url)
    x3d_available = True
except Exception:
    x3d_available = False

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not x3d_available,
                       reason='NorKyst x3d store not available from rgw.met.no'),
]

variables = ['x_sea_water_velocity', 'y_sea_water_velocity', 'x_wind', 'y_wind']


def test_norkyst_x3d_reader():
    reader = reader_netCDF_CF_generic.Reader(ds_x3d)
    assert 'stere' in reader.proj4
    for var in variables:
        assert var in reader.variables


def test_norkyst_x3d_matches_thredds():
    reader_x3d = reader_netCDF_CF_generic.Reader(ds_x3d)
    try:
        reader_thredds = reader_netCDF_CF_generic.Reader(thredds_url)
    except Exception:
        pytest.skip('NorKyst THREDDS aggregate not available')

    time = reader_x3d.start_time + timedelta(hours=48)
    lon = np.linspace(3.5, 5.5, 50)
    lat = np.linspace(59., 60.5, 50)
    kwargs = dict(variables=variables, lon=lon, lat=lat, z=np.zeros(50), time=time)
    x3d = reader_x3d.get_variables_interpolated(**kwargs)[0]
    thredds = reader_thredds.get_variables_interpolated(**kwargs)[0]

    for var in variables:
        a = np.ma.filled(np.ma.masked_invalid(x3d[var]).astype(float), np.nan)
        b = np.ma.filled(np.ma.masked_invalid(thredds[var]).astype(float), np.nan)
        assert np.isfinite(a).any()
        # int16 packing: current 1 mm/s, wind 1 cm/s
        np.testing.assert_allclose(a, b, atol=0.01, equal_nan=True)


def test_norkyst_x3d_run():
    reader = reader_netCDF_CF_generic.Reader(ds_x3d)
    o = OceanDrift(loglevel=50)
    o.add_reader(reader)
    o.seed_elements(lon=4.5, lat=59.6, radius=100, number=10,
                    time=reader.start_time + timedelta(hours=48))
    o.run(duration=timedelta(hours=6), time_step=1800)
    assert o.num_elements_active() == 10
    assert np.all(o.elements.lon != 4.5)
