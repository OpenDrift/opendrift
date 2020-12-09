import numpy as np
import pytest
import matplotlib.pyplot as plt
from tests import *
from opendrift.readers.reader_grib2 import Reader as Grib
from opendrift.readers.reader_netCDF_CF_generic import Reader as Generic
from opendrift.models.oceandrift import OceanDrift

caps_u10 = 'grb/u_vel.grb2'
caps_v10 = 'grb/v_vel.grb2'
proj4 = '+proj=stere +lat_0=90 +lat_ts=60 +lon_0=264 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +a=6371229 +b=6371229'

def test_open():
    r = Grib(caps_u10, proj4 = proj4)
    print(r)
    assert r.projected is True

def test_get_variables(benchmark):
    r = Grib(caps_u10, proj4 = proj4)

    x = [10]
    y = [60]
    z = [10]

    x, y = r.lonlat2xy(x, y)
    print(x, y)

    u = benchmark(r.get_variables, ['x_wind'], r.start_time, x, y, z)['x_wind']
    print(u)

