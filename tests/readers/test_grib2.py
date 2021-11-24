import os.path
from opendrift.readers.reader_grib2 import Reader as Grib

caps1 = 'CMC_caps_HGT_ISBL_0050_ps3km_2018090900_P005.grib2'
proj4 = '+proj=stere +lat_0=90 +lat_ts=60 +lon_0=264 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +a=6371229 +b=6371229'

def test_open(test_data):
    caps1f = os.path.join(test_data, caps1)
    r = Grib(caps1f, proj4 = proj4)
    print(r)
    assert r.projected is True

def test_get_variables(benchmark, test_data):
    caps1f = os.path.join(test_data, caps1)
    r = Grib(caps1f, proj4 = proj4)

    x = [10]
    y = [60]
    z = [10]

    x, y = r.lonlat2xy(x, y)
    print(x, y)

    u = benchmark(r.get_variables, ['geopotential_height'], r.start_time, x, y, z)['geopotential_height']
    print(u)

    assert u.shape == (2 * r.buffer , 2 * r.buffer)

