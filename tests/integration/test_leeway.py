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

import os
import time
import unittest
from datetime import datetime

from opendrift.readers import reader_basemap_landmask
from opendrift.models.leeway import Leeway


class TestLeeway(unittest.TestCase):
    """Tests for Leeway module."""

    def test_leewayprop(self):
        """Check that Leeway properties are properly read."""
        self.objectType = 85  # MED-WASTE-7
        self.lee = Leeway(loglevel=20)
        objectType = self.objectType
        self.assertEqual(self.lee.leewayprop[objectType]
                         ['Description'],
                         '>>Medical waste, syringes, small')
        self.assertEqual(self.lee.leewayprop[objectType]['DWSLOPE'], 1.79)

    def test_leeway_config_object(self):
        """Check that correct object type is fetched from config"""
        l = Leeway(loglevel=20)
        l.set_config('seed:object_type', 'Surf board with person')
        l.seed_elements(lon=4.5, lat=60, number=100,
                        time=datetime(2015, 1, 1))
        objType = l.elements_scheduled.objectType
        self.assertEqual(l.leewayprop[objType]['Description'],
                         'Surf board with person') 
        self.assertEqual(l.leewayprop[objType]['OBJKEY'],
                         'PERSON-POWERED-VESSEL-2') 

    def test_leewayrun(self):
        """Test the expected Leeway left/right split."""
        self.lee = Leeway(loglevel=30)
        self.objectType = 50  # FISHING-VESSEL-1
        self.reader_basemap = reader_basemap_landmask.Reader(
            llcrnrlon=3, llcrnrlat=59.8, projection='merc',
            urcrnrlon=6, urcrnrlat=60.5, resolution='i')
        self.lee.add_reader([self.reader_basemap])
        self.lee.seed_elements(lon=4.5, lat=60, number=100,
                               objectType=self.objectType,
                               time=datetime(2015, 1, 1))
        self.lee.fallback_values['x_wind'] = 0
        self.lee.fallback_values['y_wind'] = 10
        self.lee.fallback_values['x_sea_water_velocity'] = 0
        self.lee.fallback_values['y_sea_water_velocity'] = 0
        # Check that 7 out of 100 elements strand towards coast
        self.lee.run(steps=24, time_step=3600)
        self.assertEqual(self.lee.num_elements_scheduled(), 0)
        self.assertEqual(self.lee.num_elements_active(), 97)
        self.assertEqual(self.lee.num_elements_deactivated(), 3)  # stranded
        self.lee.export_ascii('leeway_ascii.txt')
        os.remove('leeway_ascii.txt')


if __name__ == '__main__':
    unittest.main()
