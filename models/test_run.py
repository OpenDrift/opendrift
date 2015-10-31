#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from datetime import datetime, timedelta
import os

import numpy as np

from readers import reader_ArtificialOceanEddy
from readers import reader_basemap_landmask
from models.oceandrift import OceanDrift


class TestArray(unittest.TestCase):
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
        self.o.seed_point(2.0, 61.0, radius=0, number=9,
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

    def test2_seed_points(self):
        """Test a model run"""
        self.make_OceanDrift_object()
        self.o.seed_point([2.0, 4.5, 3.0], [61.0, 60.0, 62.0],
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

if __name__ == '__main__':
    unittest.main()
