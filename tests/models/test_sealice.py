import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.sealice import SeaLice

def test_sealice():
    o = SeaLice(loglevel=30)
    reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF_20160202_subset.nc')
    o.add_reader([reader_arome])
    lat = 67.711251; lon = 13.556971  # Lofoten
    o.seed_elements(lon, lat, radius=5000, number=1000,
                    time=reader_arome.start_time)
    o.run(steps=24, time_step=3600)
    np.testing.assert_almost_equal(o.elements.lon.max(), 15.864, 2)
