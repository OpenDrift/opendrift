import os.path
from netCDF4 import Dataset

def test_netcdf_segfault(tmpdir):
    # https://github.com/Unidata/netcdf4-python/issues/982
    d1 = Dataset(os.path.join(tmpdir, 'test.nc'), 'w')
    d1.close()
    d2 = Dataset('https://thredds.met.no/thredds/dodsC/meps25epsarchive/2019/11/26/meps_mbr0_extracted_2_5km_20191126T00Z.nc')
    print(d2.variables['time'])

