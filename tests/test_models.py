#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import unittest
from datetime import datetime

from opendrift.readers import reader_ArtificialOceanEddy
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil
from opendrift.models.windblow import WindBlow


class TestArray(unittest.TestCase):
    """Tests for OpenDrift models"""

    def setUp(self):
        self.o = OpenOil()
        self.fake_eddy = reader_ArtificialOceanEddy.Reader(2, 62)
        self.fake_eddy.start_time = datetime(2015, 1, 1)
        self.o.add_reader([self.fake_eddy])

    def test_seed(self):
        """Test seeding"""
        self.o.seed_elements(lon=4, lat=60, number=100,
                             time=self.fake_eddy.start_time)
        self.assertEqual(len(self.o.elements), 0)
        self.assertEqual(len(self.o.elements_deactivated), 0)
        self.assertEqual(len(self.o.elements_scheduled), 100)

    def test_windblow(self):
        o = WindBlow(loglevel=20)
        reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF.nc_20160202_subset')
        o.add_reader([reader_arome])
        lat = 67.711251; lon = 13.556971  # Lofoten
        o.seed_elements(lon, lat, radius=5000, number=1000,
                        time=reader_arome.start_time)
        o.run(steps=48*4, time_step=900)
        self.assertAlmostEqual(o.elements.lon.max(), 17.64, 2)


if __name__ == '__main__':
    unittest.main()
