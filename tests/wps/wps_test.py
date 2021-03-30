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
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import unittest
import os
from datetime import datetime, timedelta

import numpy as np

from opendrift.models.oceandrift import OceanDrift
from opendrift.models.leeway import Leeway
from opendrift.models.openoil import OpenOil
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_constant
from opendrift.readers import reader_from_url



class TestWPS(unittest.TestCase):
    """Tests for wps simulations"""

    def test_leeway_today(self):
        o = Leeway(loglevel=0)
        o.add_readers_from_file(o.test_data_folder() +
            '../../opendrift/scripts/data_sources.txt')
        o.seed_elements(lon=14, lat=67.85, number=100, radius=1000,
                        time=datetime.now())
        o.run(steps=15)
        print (o)
        self.assertEqual(o.steps_calculation, 15)

    def test_leeway_yesterday(self):
        o = Leeway(loglevel=0)
        o.add_readers_from_file(o.test_data_folder() +
            '../../opendrift/scripts/data_sources.txt')
        o.seed_elements(lon=14, lat=67.85, number=100, radius=1000,
                        time=datetime.now() - timedelta(days=1))
        o.run(steps=15)
        o.export_ascii('leeway_ascii.txt')
        os.remove('leeway_ascii.txt')
        print (o)
        self.assertEqual(o.steps_calculation, 15)

    def test_leeway_global_today(self):
        o = Leeway(loglevel=0)
        o.add_readers_from_file(o.test_data_folder() +
            '../../opendrift/scripts/data_sources.txt')
        o.seed_elements(lon=50, lat=29, number=100, radius=1000,
                        time=datetime.now())
        o.run(steps=15)
        print (o)
        self.assertEqual(o.steps_calculation, 15)

    def test_leeway_global_one_month_ago(self):
        o = Leeway(loglevel=0)
        o.add_readers_from_file(o.test_data_folder() +
            '../../opendrift/scripts/data_sources.txt')
        o.seed_elements(lon=50, lat=29, number=100, radius=1000,
                        time=datetime.now() - timedelta(days=30))
        o.run(steps=15)
        o.export_ascii('leeway_ascii.txt')
        os.remove('leeway_ascii.txt')
        print (o)
        self.assertEqual(o.steps_calculation, 15)

    def test_openoil_today(self):
        o = OpenOil(loglevel=0)
        o.add_readers_from_file(o.test_data_folder() +
            '../../opendrift/scripts/data_sources.txt')
        o.seed_elements(lon=14, lat=67.85, number=100, radius=1000,
                        time=datetime.now())
        o.run(steps=15)
        print (o)
        self.assertEqual(o.steps_calculation, 15)

    def test_openoil_yesterday(self):
        o = OpenOil(loglevel=0)
        o.add_readers_from_file(o.test_data_folder() +
            '../../opendrift/scripts/data_sources.txt')
        o.seed_elements(lon=14, lat=67.85, number=100, radius=1000,
                        time=datetime.now() - timedelta(days=1))
        o.run(steps=15)
        print (o)
        self.assertEqual(o.steps_calculation, 15)

    def test_openoil_global_today(self):
        o = OpenOil(loglevel=0)
        o.add_readers_from_file(o.test_data_folder() +
            '../../opendrift/scripts/data_sources.txt')
        o.seed_elements(lon=50, lat=29, number=100, radius=1000,
                        time=datetime.now())
        o.run(steps=15)
        print (o)
        self.assertEqual(o.steps_calculation, 15)

    def test_openoil_global_one_month_ago(self):
        o = OpenOil(loglevel=0)
        o.add_readers_from_file(o.test_data_folder() +
            '../../opendrift/scripts/data_sources.txt')
        o.seed_elements(lon=50, lat=29, number=100, radius=1000,
                        time=datetime.now() - timedelta(days=30))
        o.run(steps=15)
        print (o)
        self.assertEqual(o.steps_calculation, 15)


    #def test_oildrift_backwards(self):
    #    o = OpenOil(loglevel=20)
    #    reader_constant_wind = \
    #        reader_constant.Reader({'x_wind':5, 'y_wind': 6})
    #    o.add_reader(reader_constant_wind)

    #    o.add_readers_from_list(reader_list, lazy=True)

    #    self.assertEqual(len(o._lazy_readers()), 4)
    #    o.seed_elements(lon=14, lat=67.85,
    #                    time=datetime(2016, 2, 2, 12))
    #    o.set_config()
    #    o.run(steps=5)
    #    self.assertEqual(len(o._lazy_readers()), 2)
    #    self.assertEqual(len(o.discarded_readers), 1)

    #def test_lazy_reader_oildrift_real(self):
    #    o = OpenOil(loglevel=0)
    #    o.add_readers_from_file(o.test_data_folder() +
    #        '../../opendrift/scripts/data_sources.txt')

    #    o.seed_elements(lon=4, lat=60.0,
    #                    time=datetime(2018, 7, 2, 12))
    #    o.run(steps=5)
    #    print (o)

if __name__ == '__main__':
    unittest.main()
