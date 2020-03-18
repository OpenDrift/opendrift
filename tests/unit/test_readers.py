import unittest
from datetime import datetime, timedelta
import numpy as np
import pytest
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_timeseries
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

    def test_reader_timeseries(self):
        r = reader_timeseries.Reader({
            'time': [datetime(2020, 1, 1, 0),
                     datetime(2020, 1, 1, 6),
                     datetime(2020, 1, 1, 12),
                     datetime(2020, 1, 1, 18),
                     datetime(2020, 1, 2, 0)],
            'x_wind': [3, 3, 5, 8, 6],
            'y_wind': [1, 2, 4, 3, 2]})
        d = r.get_variables(['x_wind', 'y_wind'],
                            x = np.array([2, 3, 4]),
                            y = np.array([58, 59, 60]),
                            time=datetime(2020, 1, 1, 12))
        self.assertEqual(d['x_wind'][1], 5)
        d = r.get_variables(['x_wind', 'y_wind'],
                            x = np.array([2, 3, 4]),
                            y = np.array([58, 59, 60]),
                            time=datetime(2020, 1, 1, 17))
        self.assertEqual(d['x_wind'][1], 8)
        # Insert invalid data
        self.assertRaises(ValueError,
            reader_timeseries.Reader, {
                'time': [datetime(2020, 1, 1, 0),
                         datetime(2020, 1, 1, 6),
                         datetime(2020, 1, 1, 12),
                         datetime(2020, 1, 1, 18),
                         datetime(2020, 1, 2, 0)],
                'x_wind': [3, 3],
                'y_wind': [1, 2, 4, 3, 2]})

if __name__ == '__main__':
    unittest.main()
