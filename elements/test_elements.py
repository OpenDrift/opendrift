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

import numpy as np

from elements import LagrangianArray


#class TestScalar(unittest.TestCase):
#    """Tests for LagrangianArray scalar"""
#
#    def test_init(self):
#        """Initiate a Lagrangian scalar, with default"""
#        e = LagrangianArray(lon=4.0, lat=60.0)
#        self.assertEqual(e.shape, ())
#        self.assertRaises(TypeError, len, e)  # scalar has no length
#
#    def test_allargs(self):
#        """Initiate a Lagrangian scalar, with all arguments"""
#        e = LagrangianArray(lon=4.0, lat=60.0, depth=10.0, active=False)
#        self.assertEqual(e.lon, 4.0)
#        self.assertEqual(e.lat, 60.0)
#        self.assertEqual(e.depth, 10.0)
#        self.assertEqual(e.active, False)
#
#    def test_defaults(self):
#        """Check that the defaults are working properly"""
#        e = LagrangianArray(lon=4.0, lat=60.0)
#        self.assertEqual(e.depth, 0.0)
#        self.assertEqual(e.active, True)
#
#    def test_no_index(self):
#        """Elements can't be indexed, e[0] raises IndexError"""
#        e = LagrangianArray(lon=4.0, lat=60.0)
#        self.assertRaises(IndexError, e.__getitem__, 0)
#
#    # Not implemented
#    def rest_tuple(self):
#        """Should be able to take the tuple"""
#        e = LagrangianArray(lon=4.0, lat=60.0)
#        self.assertEqual(tuple(e), (4.0, 60.0, 0.0, True))
#        self.assertEqual(LagrangianArray(*tuple(e)), e)
#
#    # Not implemented
#    def rest_str(self):
#        e = LagrangianArray(lon=4.0, lat=60.0, depth=10.0)
#        self.assertEqual(str(e), "(4.0, 60.0, 10.0, True)")


class TestArray(unittest.TestCase):
    """Tests for (non-scalar) LagrangianArray"""

    def test_init(self):
        """Chcek that init works"""
        A = LagrangianArray(lon=4.0, lat=60.0)
        self.assertEqual(len(A), 1)

    def test_fromarrays(self):
        """Make LagrangianArray from attribute arrays"""
        lon = [4.0, 4.5, 5.0]
        lat = 60.0
        depth = 10
        A = LagrangianArray(lon=lon, lat=lat, depth=depth)
        self.assertEqual(len(A), 3)
        self.assertEqual(list(A.lon), [4.0, 4.5, 5.0])

    def test_empty(self):
        """Empty LagrangianArray works as expected"""
        A0 = LagrangianArray()
        self.assertEqual(len(A0), 0)
        self.assertFalse(A0)
        self.assertFalse(A0.lon)

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

#    def test_indexing(self):
#        """LagrangianArray items are Lagrangian elements"""
#        e1 = LagrangianArray(lon=2, lat=60, depth=10)
#        e2 = LagrangianArray(lon=2, lat=61, depth=10)
#        e3 = LagrangianArray(lon=3, lat=61)
#        A = LagrangianArray.fromelements([e1, e2, e3])
#        a = A[1]
#        self.assertIsInstance(a, LagrangianArray)
#        self.assertEqual(a.shape, ())   # a is scalar
#        self.assertEqual(a.lat, e2.lat)
#        self.assertEqual(a.lat, A.lat[1])
#        # A[len(A)] raises IndexError
#        self.assertRaises(IndexError, A.__getitem__, len(A))
#
#    def test_iteration(self):
#        """A LagrangianArray is an iterable"""
#        e1 = LagrangianArray(lon=2, lat=60, depth=10)
#        e2 = LagrangianArray(lon=2, lat=61)
#        e3 = LagrangianArray(lon=3, lat=61)
#        A = LagrangianArray.fromelements([e1, e2, e3])
#        for i, e in enumerate(A):
#            self.assertEqual(A[i].lon, e.lon)
#
#    def test_setitem(self):
#        """Setting items can mutate a LagrangianArray"""
#        e = LagrangianArray(lon=2, lat=66)
#        lon = [4.2, 4.4, 4.6]
#        lat = 60.0
#        depth = 100.0
#        A = LagrangianArray.fromarrays(lon=lon, lat=lat, depth=depth)
#        idA = id(A)
#        A[1] = e                      # setitem
#        self.assertEqual(A.lat[1], 66)  # A has changed
#        self.assertEqual(id(A), idA)  # A is still A
#
#    def test_append(self):
#        """Elements can be appended to a LagrangianArray"""
#        e1 = LagrangianArray(lon=2, lat=60, depth=10)
#        e2 = LagrangianArray(lon=2, lat=61)
#        e3 = LagrangianArray(lon=3, lat=61)
#        A = LagrangianArray.fromelements([e1, e2, e3])
#        idA = id(A)
#        lenA = len(A)
#        shapeA = A.shape
#        e = LagrangianArray(lon=2.0, lat=6.0, depth=100)  # station M
#        A.append(e)
#        # The length should increase by one
#        self.assertEqual(len(A), lenA+1)
#        # The shape should be adjusted
#        self.assertEqual(A.shape[0], shapeA[0] + 1)
#        #  The last element should be e
#        self.assertEqual(A[-1].lon, e.lon)
#        # A is still A
#        self.assertEqual(id(A), idA)
#
#    def test_extend(self):
#        """Test concatenation"""
#        e1 = LagrangianArray(lon=2, lat=60)
#        e2 = LagrangianArray(lon=2, lat=61)
#        e3 = LagrangianArray(lon=3, lat=61)
#        A = LagrangianArray.fromelements([e1, e2, e3])
#        idA = id(A)
#        lenA = len(A)
#        shapeA = A.shape
#        lon = [4, 5]
#        lat = [59, 58]
#        B = LagrangianArray.fromarrays(lon=lon, lat=lat)
#        A.extend(B)
#        # Length OK
#        self.assertEqual(len(A), lenA + len(B))
#        # Shape OK
#        self.assertEqual(A.shape[0], shapeA[0] + B.shape[0])
#        # First new element is correct
#        self.assertEqual(A[lenA].lon, B[0].lon)
#        # Last element is correct
#        self.assertEqual(A[-1].lat, B[-1].lat)
#        # A is still A
#        self.assertEqual(id(A), idA)

if __name__ == '__main__':
    unittest.main()
