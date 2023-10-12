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
# Copyright 2015,2020 Knut-Frode Dagestad, MET Norway

import unittest
import pytest
from datetime import datetime, timedelta
import os
import numpy as np
import opendrift
from opendrift.readers import reader_oscillating
from opendrift.models.oceandrift import OceanDrift


class TestXarray(unittest.TestCase):
    """Tests for Xarray / dask functionality"""

    def test_density(self):
        """Test density"""
        outfile = 'test_xarray.nc'
        if os.path.exists(outfile):
            os.remove(outfile)
        o = OceanDrift(loglevel=20)
        t1 = datetime.now()
        t2 = t1 + timedelta(hours=6)
        reader_x = reader_oscillating.Reader('x_sea_water_velocity',
                        amplitude=1, zero_time=t1)
        reader_y = reader_oscillating.Reader('y_sea_water_velocity',
                        amplitude=1, zero_time=t2)
        o.add_reader([reader_x, reader_y])
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('drift:horizontal_diffusivity', 10)
        o.seed_elements(time=t1, lon=4, lat=60, number=100,
                        origin_marker=0)
        o.seed_elements(time=[t1, t2], lon=4.2, lat=60.2, number=100,
                        origin_marker=1)
        o.seed_elements(time=[t1, t2], lon=4.1, lat=60.1, number=100,
                        origin_marker=2)
        o.run(duration=timedelta(hours=12), time_step=1800, outfile=outfile)
        #o.plot(fast=True)
        density_pixelsize_m=5000
        H, Hsub, Hsurf, lon_array, lat_array = o.get_density_array(pixelsize_m=density_pixelsize_m)

        ox = opendrift.open_xarray(outfile)
        Hx = ox.get_histogram(pixelsize_m=density_pixelsize_m)
        self.assertAlmostEqual(lon_array[0], 3.94, 1)
        self.assertAlmostEqual(lon_array[-1], 4.76, 1)
        self.assertAlmostEqual(Hx.lon_bin[0].values, 3.90, 1)
        self.assertAlmostEqual(Hx.lon_bin[-1].values, 4.67, 1)
        self.assertEqual(Hx.sum(dim='origin_marker').shape, H.shape)
        Hsum = H.sum(axis=1).sum(axis=1)
        Hxsum = Hx.sum(('lon_bin', 'lat_bin', 'origin_marker'))
        self.assertEqual(Hsum[0], 118)
        self.assertEqual(Hxsum[0], 118)
        self.assertEqual(Hsum[-1], 300)
        self.assertEqual(Hxsum[-1], 300)
        os.remove(outfile)

if __name__ == '__main__':
    unittest.main()
