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

import numpy as np

from opendrift.elements.elements import LagrangianArray


class TestElements(unittest.TestCase):
    """Tests for (non-scalar) LagrangianArray"""

    def test_init(self):
        """Chcek that init works"""
        A = LagrangianArray(lon=4.0, lat=60.0)
        self.assertEqual(len(A), 1)

    def test_fromarrays(self):
        """Make LagrangianArray from attribute arrays"""
        lon = [4.0, 4.5, 5.0]
        lat = 60.0
        z = -10
        A = LagrangianArray(lon=lon, lat=lat, z=z)
        self.assertEqual(len(A), 3)
        self.assertEqual(list(A.lon), [4.0, 4.5, 5.0])

    def test_empty(self):
        """Empty LagrangianArray works as expected"""
        A0 = LagrangianArray()
        self.assertEqual(len(A0), 0)
        self.assertFalse(A0)
        self.assertTrue(len(A0.lon)==0)

    def test_move(self):
        A1 = LagrangianArray(lon = [4.1], lat = [60.1])
        A2 = LagrangianArray(lon = [4.0, 4.5], lat = [60.0, 60.5])
        A3 = LagrangianArray(lon = [1.0, 1.5, 1.7], lat = [60.0, 60.5, 61.0])
        # Moving first and third element from A3 to A2
        A3.move_elements(A2, np.array([True, False, True], dtype=bool))
        self.assertEqual(len(A3), 1)
        self.assertEqual(len(A2), 4)
        A2.move_elements(A3, np.array([False, True, False, False], dtype=bool))
        self.assertEqual(len(A2), 3)
        self.assertEqual(len(A1), 1)
        A1.move_elements(A3, np.array([True], dtype=bool))
        self.assertEqual(len(A1), 0)

    def test_extend(self):
        """Test concatenation"""
        e1 = LagrangianArray(lon=2, lat=60, z=[0, 1, 2])
        e2 = LagrangianArray(lon=2, lat=61, z=-1)
        e3 = LagrangianArray(lon=3, lat=61, z=0)

        self.assertTrue(hasattr(e1.z, '__len__'))
        self.assertFalse(hasattr(e2.z, '__len__'))
        # Extend 
        e1.extend(e2)
        e2.extend(e3)
        # Check that identical scalars remain scalars, arrays are concatenated
        self.assertTrue(hasattr(e1.z, '__len__'))
        self.assertTrue(hasattr(e2.z, '__len__'))
        self.assertListEqual(list(e1.z), [0, 1, 2, -1])
        self.assertListEqual(list(e2.lon), [2., 3.])
        self.assertEqual(e2.lat, 61.0)

        self.assertEqual(len(e1), 4)
        self.assertEqual(len(e2), 2)
        self.assertEqual(len(e3), 1)

if __name__ == '__main__':
    unittest.main()
