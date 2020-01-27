import unittest
from datetime import timedelta
import numpy as np
import pytest
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

class TestReaders(unittest.TestCase):

    def test_timeseries_at_position(self):

        o = OceanDrift()
        reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF.nc_20160202_subset')

        ts = reader_arome.get_timeseries_at_position(
            lon=10, lat=68, variables = ['x_wind', 'y_wind'])

        assert len(ts['time']) == 49
        x_wind = ts['x_wind']
        self.assertAlmostEqual(x_wind[0], -1.081, 2)
        self.assertAlmostEqual(x_wind[-1], -3.991, 2)

if __name__ == '__main__':
    unittest.main()
