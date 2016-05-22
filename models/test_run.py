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
import os
import inspect

import numpy as np

from readers import reader_ArtificialOceanEddy
from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.oceandrift import OceanDrift

script_folder = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))

class TestRun(unittest.TestCase):
    """Tests for (non-scalar) LagrangianArray"""

    def make_OceanDrift_object(self):
        self.o = OceanDrift(loglevel=20)
        self.fake_eddy = reader_ArtificialOceanEddy.Reader(2, 62)
        self.o.use_block = False
        self.o.config['runge_kutta'] = False
        self.reader_basemap = reader_basemap_landmask.Reader(
            llcrnrlon=-1.5, llcrnrlat=59,
            urcrnrlon=7, urcrnrlat=64, resolution='i')
        self.o.add_reader([self.fake_eddy, self.reader_basemap])

    def test1_seed_single_point_over_time(self):
        """Test a model run"""
        self.make_OceanDrift_object()
        self.o.seed_elements(2.0, 61.0, radius=0, number=9,
                             time=[datetime(2015, 1, 1), datetime(2015, 1, 3)])

        # Check that 6 elements are scheduled, but none seeded
        self.assertEqual(self.o.num_elements_scheduled(), 9)
        self.assertEqual(self.o.num_elements_active(), 0)
        self.assertEqual(self.o.num_elements_activated(), 0)
        self.assertEqual(self.o.num_elements_deactivated(), 0)
        self.assertEqual(self.o.num_elements_total(), 9)
        # Run simulation
        self.o.run(steps=30, outfile='unittest.nc')
        # Check that 1 element is deactivated (stranded),
        # 1 yet not seeded and 4 active
        self.assertEqual(self.o.num_elements_scheduled(), 4)
        self.assertEqual(self.o.num_elements_active(), 5)
        self.assertEqual(self.o.num_elements_activated(), 5)
        self.assertEqual(self.o.num_elements_deactivated(), 0)
        self.assertEqual(self.o.num_elements_total(), 9)

    def test2_seed_elementss(self):
        """Test a model run"""
        self.make_OceanDrift_object()
        self.o.seed_elements([2.0, 4.5, 3.0], [61.0, 60.0, 62.0],
                             radius=0, number=9,
                             time=[datetime(2015, 1, 1), datetime(2015, 1, 3)])

        # Check that 6 elements are scheduled, but none seeded
        self.assertEqual(self.o.num_elements_scheduled(), 9)
        self.assertEqual(self.o.num_elements_active(), 0)
        self.assertEqual(self.o.num_elements_activated(), 0)
        self.assertEqual(self.o.num_elements_deactivated(), 0)
        self.assertEqual(self.o.num_elements_total(), 9)
        # Run simulation
        self.o.run(steps=30, outfile='unittest.nc')
        # Check that 1 element is deactivated (stranded),
        # 1 yet not seeded and 4 active
        self.assertEqual(self.o.num_elements_scheduled(), 4)
        self.assertEqual(self.o.num_elements_active(), 4)
        self.assertEqual(self.o.num_elements_activated(), 5)
        self.assertEqual(self.o.num_elements_deactivated(), 1)
        self.assertEqual(self.o.num_elements_total(), 9)

    def test3_run_import(self):
        """Import output file from previous test, and check elements"""
        self.o = OceanDrift(loglevel=20)
        self.o.io_import_file('unittest.nc')
        self.assertEqual(self.o.num_elements_active(), 4)
        self.assertEqual(self.o.num_elements_activated(), 5)
        self.assertEqual(self.o.num_elements_deactivated(), 1)
        self.assertEqual(self.o.num_elements_total(), 5)

    def test4_cleaning(self):
        """Cleaning up"""
        os.remove('unittest.nc')

    def test_output_time_step(self):
        o1 = OceanDrift(loglevel=20)
        norkyst = reader_netCDF_CF_generic.Reader(script_folder + '/../test_data/16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        basemap = reader_basemap_landmask.Reader(
            llcrnrlon=4, llcrnrlat=59.8, urcrnrlon=6, urcrnrlat=61,
            resolution='h', projection='merc')
        o1.add_reader([basemap, norkyst])
        o1.seed_elements(5.25, 60.2, radius=3000, number=100,
                        time=norkyst.start_time)
        o1.run(duration=timedelta(hours=12),
                   time_step=timedelta(minutes=30),
                   time_step_output=timedelta(minutes=30),
                   outfile='test_time_step30.nc')
        # Check length of time array and output array
        time = o1.get_time_array()[0]
        self.assertEqual(o1.history.shape[1], len(time))
        self.assertEqual(o1.start_time, time[0])
        self.assertEqual(o1.time, time[-1])
        # Second run, with larger output time step
        o2 = OceanDrift(loglevel=20)
        o2.add_reader([basemap, norkyst])
        o2.seed_elements(5.25, 60.2, radius=3000, number=100,
                        time=norkyst.start_time)
        o2.run(duration=timedelta(hours=12),
                   time_step=timedelta(minutes=30),
                   time_step_output=timedelta(minutes=60),
                   outfile='test_time_step60.nc')
        self.assertEqual(o1.history.shape, (100,25))
        self.assertEqual(o2.history.shape, (100,13))
        # Check that start and end conditions (longitudes) are idential
        self.assertItemsEqual(o1.history['lon'][:,24].compressed(), 
                              o2.history['lon'][:,12].compressed())
        self.assertItemsEqual(o1.history['lon'][:,0].compressed(),
                              o2.history['lon'][:,0].compressed())
        # Check that also run imported from file is identical
        o1i = OceanDrift(loglevel=20)
        o1i.io_import_file('test_time_step30.nc')
        o2i = OceanDrift(loglevel=20)
        o2i.io_import_file('test_time_step60.nc')
        os.remove('test_time_step30.nc')
        os.remove('test_time_step60.nc')
        self.assertItemsEqual(o2i.history['lon'][:,12].compressed(),
                              o2.history['lon'][:,12].compressed())
        # Check number of activated elements
        self.assertEqual(o1.num_elements_total(), o2.num_elements_total())
        self.assertEqual(o1.num_elements_total(), o1i.num_elements_total())
        self.assertEqual(o1.num_elements_total(), o2i.num_elements_total())
        # Check number of deactivated elements
        self.assertEqual(o1.num_elements_deactivated(),
                         o2.num_elements_deactivated())
        self.assertEqual(o1.num_elements_deactivated(),
                         o1i.num_elements_deactivated())
        self.assertEqual(o1.num_elements_deactivated(),
                         o2i.num_elements_deactivated())


if __name__ == '__main__':
    unittest.main()
