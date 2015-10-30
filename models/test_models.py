#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
from datetime import datetime, timedelta

import numpy as np

from readers import reader_ArtificialOceanEddy
from models.openoil import OpenOil


class TestArray(unittest.TestCase):
    """Tests for (non-scalar) LagrangianArray"""

    def setUp(self):
        self.o = OpenOil()
        self.fake_eddy = reader_ArtificialOceanEddy.Reader(2, 62)
        self.fake_eddy.start_time = datetime(2015, 1, 1)
        self.o.add_reader([self.fake_eddy])

    def test_seed(self):
        """Test seeding"""
        self.o.seed_point(lon=4, lat=60, number=100,
                          time=self.fake_eddy.start_time)
        self.assertEqual(len(self.o.elements), 0)
        self.assertEqual(len(self.o.elements_deactivated), 0)
        self.assertEqual(len(self.o.elements_scheduled), 100)



if __name__ == '__main__':
    unittest.main()
