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
        self.o.seed_elements(lon=4, lat=60, number=100,
                             time=self.fake_eddy.start_time)
        self.assertEqual(len(self.o.elements), 0)
        self.assertEqual(len(self.o.elements_deactivated), 0)
        self.assertEqual(len(self.o.elements_scheduled), 100)


if __name__ == '__main__':
    unittest.main()
