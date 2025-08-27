import pytest
import numpy as np
from datetime import datetime, timedelta
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift

def test_croco(test_data):

    reader_croco = reader_ROMS_native.Reader(test_data + 'croco/croco_his.nc',
        gridfile = test_data + 'croco/croco_grd.nc')

    o = OceanDrift(loglevel=50)
    o.add_reader(reader_croco)
    o.seed_elements(lon=10, lat=-30, time=reader_croco.start_time)
    o.run(end_time=reader_croco.end_time, time_step=timedelta(hours=3))
    np.testing.assert_almost_equal(o.elements.lon, 9.99359, 5)
    np.testing.assert_almost_equal(o.elements.lat, -30.0085, 5)
