import numpy as np
from opendrift.readers import reader_schism_native

proj4str_nztm = '+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'


def test_schicm_load():
    schism_native = reader_schism_native.Reader(
        filename=
        'https://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/netcdf_unstructured_samples/schism_marl20080101_00z_3D.nc',
        proj4=proj4str_nztm,
        use_3d=True)

    print(schism_native)

def test_get_variables():
    schism_native = reader_schism_native.Reader(
        filename=
        'https://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/netcdf_unstructured_samples/schism_marl20080101_00z_3D.nc',
        proj4=proj4str_nztm,
        use_3d=True)
    e, _ = schism_native.get_variables_interpolated(['x_sea_water_velocity'], time = schism_native.start_time, lon = 173., lat = -40.)
    print(e)

    np.testing.assert_almost_equal(e['x_sea_water_velocity'], [-0.16828385])

