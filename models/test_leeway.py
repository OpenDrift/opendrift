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
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt

from readers import reader_basemap_landmask
from models.leeway import Leeway
from models.windblow import WindBlow


class TestArray(unittest.TestCase):
    """Tests for Leeway module."""

    def setUp(self):
        self.objectType = 50  # FISHING-VESSEL-1
        self.lee = Leeway(loglevel=20)
        #print self.lee.leewayprop.values()[0]

        #self.lee = WindBlow(loglevel=0)
        self.reader_basemap = reader_basemap_landmask.Reader(
            llcrnrlon=3, llcrnrlat=59, projection='merc',
            urcrnrlon=6, urcrnrlat=61, resolution='i')
        self.lee.add_reader([self.reader_basemap])
        self.lee.fallback_values['x_wind'] = 0
        self.lee.fallback_values['y_wind'] = 10
        self.lee.fallback_values['x_sea_water_velocity'] = 0
        self.lee.fallback_values['y_sea_water_velocity'] = 0

    def test_leewayprop(self):
        """Check that Leeway properties are properly read."""
        objectType = self.objectType
        self.assertEqual(self.lee.leewayprop[objectType]
            ['Description'], ' Fishing vessel, general (mean values)\n')
        self.assertEqual(self.lee.leewayprop[objectType]['DWSLOPE'], 2.47)

    def test_leewayrun(self):
        """Test the expected Leeway left/right split."""
        self.lee.seed_elements(lon=4.5, lat=60, number=100,
                             objectType=self.objectType,
                             time=datetime(2015, 1, 1))
        # Check that 7 out of 100 elements strand towards coast
        self.lee.run(steps=24, time_step=3600, outfile='leewaytest.nc')
        self.assertEqual(self.lee.num_elements_scheduled(), 0)
        self.assertEqual(self.lee.num_elements_active(), 98)
        self.assertEqual(self.lee.num_elements_deactivated(), 2)  # stranded


if __name__ == '__main__':
    unittest.main()
