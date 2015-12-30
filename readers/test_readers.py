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
import os

import numpy as np

import reader_netCDF_CF_generic
import reader_ROMS_native


class TestReaders(unittest.TestCase):
    """Tests for readers"""

    def test_reader_netcdf(self):
        """Check reader functionality."""
        readers = [reader_netCDF_CF_generic.Reader(
            'test_data/norkyst800_subset_16Nov2015.nc')]
        if os.path.exists('/disk2/data/roms/ocean_his_0170.nc'):
            readers.append(reader_ROMS_native.Reader(
                '/disk2/data/roms/ocean_his_0170.nc'))
        for r in readers:
            print r
            # Make four points:
            #  1) outside lower left, 2) lower left,  3) center of domain
            #  4) outside upper right
            # and assure that only 2) and 3) are marked as covered
            # Upper right is skipped, as lonlat2xy may lie slightly outside
            x = np.array([r.xmin - r.delta_x, r.xmin, (r.xmin + r.xmax)/2,
                          r.xmax + r.delta_x])
            y = np.array([r.ymin - r.delta_y, r.ymin, (r.ymin + r.ymax)/2,
                          r.ymax + r.delta_y])
            lons, lats = r.xy2lonlat(x,  y)
            covered = r.covers_positions(lons, lats, 0)
            self.assertEqual(covered.tolist(), [1, 2])

            self.assertTrue(r.covers_time(r.start_time))
            self.assertFalse(r.covers_time(r.start_time - r.time_step))
            self.assertFalse(r.proj.is_latlong())


if __name__ == '__main__':
    unittest.main()
