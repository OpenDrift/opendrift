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
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

import reader_netCDF_CF_generic
import reader_ROMS_native
from interpolation import ReaderBlock

class TestInterpolation(unittest.TestCase):
    """Tests spatial interpolation"""


    def get_synthetic_data_dict(self):
        data_dict = {}
        data_dict['x'] = np.linspace(-70, 470, 200)
        data_dict['y'] = np.linspace(10, 340, 100)
        data_dict['z'] = np.array([-0, -3, -10, -25, -100])
        # Make a horizontal slice
        xg, yg = np.meshgrid(data_dict['x'], data_dict['y'])
        slice1 = np.ma.array(np.cos(np.radians(xg)) +
                             np.sin(np.radians(yg)))
        # Add some holes
        slice1[0:40, 50:60] = np.nan
        slice1[40:60, 100:120] = np.nan
        slice1[20:22, 30:32] = np.nan
        slice1 = np.ma.masked_invalid(slice1)

        # Make another horizontal slice ("below") with more holes
        slice2 = slice1*1.1
        slice2[70:80, 20:28] = np.nan
        
        # Add a 2D and a 3D variable to dictionary
        data_dict['var2d'] = slice1
        data_dict['var3d'] = np.ma.array([slice1, slice2, 1.2*slice1,
                                          1*3*slice1, 10*slice1])
        data_dict['time'] = datetime.now()

        # Generate some points
        x = np.linspace(data_dict['x'].min(), data_dict['x'].max(), 100)
        y = np.linspace(data_dict['y'].min(), data_dict['y'].max(), 100)

        return data_dict, x, y

    def test_interpolation_synthetic_data(self):

        data_dict, x, y = self.get_synthetic_data_dict()
        # Make block from dictionary, and apply tests
        b = ReaderBlock(data_dict, interpolation_method='ndimage')
        self.assertEqual(b.data_dict['var2d'].shape,
                         (len(b.y), len(b.x)))
        self.assertEqual(b.data_dict['var3d'].shape,
                         (len(b.z), len(b.y), len(b.x)))
        # Make some element positions
        #plt.pcolor(b.x, b.y, data_dict['var2d'])
        #plt.scatter(x, y)
        #plt.show()
        preparation = b.interpolate(x, y, variables='prepare')
        values = b.interpolate_horizontal_layers(data_dict['var2d'])
        # Check that 15 values are masked when using ndimage
        self.assertEqual(sum(values.mask), 15)
        # 3D interpolation
        values = b.interpolate_horizontal_layers(data_dict['var3d'])
        self.assertEqual(values.shape, (len(b.z), len(x)))
        self.assertEqual(sum(values.ravel().mask), 75)

        b = ReaderBlock(data_dict, interpolation_method='linearND')
        preparation = b.interpolate(x, y, variables='prepare')
        # Check that linearND interpolate into 'holes'
        values = b.interpolate_horizontal_layers(data_dict['var2d'])
        self.assertEqual(sum(values.mask), 0)
        values = b.interpolate_horizontal_layers(data_dict['var3d'])
        self.assertEqual(sum(values.ravel().mask), 0)

    def test_interpolation_synthetic_data_main(self):

        data_dict, x, y = self.get_synthetic_data_dict()
        # Make block from dictionary, and apply tests
        b = ReaderBlock(data_dict, interpolation_method='ndimage')
        env = b.interpolate(x, y)
        self.assertEqual(env['var2d'].shape, (len(x),))
        self.assertEqual(env['var3d'].shape, (len(data_dict['z']), len(x)))

        # Check values of a vertical transect
        self.assertTrue(np.allclose(env['var3d'][:,10],
            [1.64879798585, 1.81367778444, 1.97855758302,
             4.94639395756, 16.4879798585]))

    def atest_interpolation(self):
        """Test interpolation."""
        reader = reader_netCDF_CF_generic.Reader(
            '../test_data/norkyst800_subset_16Nov2015.nc')

        print reader
        num_points = 1000
        x = np.random.uniform(reader.xmin, reader.xmax, num_points)
        y = np.random.uniform(reader.ymin, reader.ymax, num_points)
        lon, lat = reader.xy2lonlat(x, y)
        variables = ['x_sea_water_velocity', 'y_sea_water_velocity']

        # Read a block of data
        data = reader.get_variables(variables, time=reader.start_time,
                                     x=x, y=y, z=0, block=True)
        slice2d = data[variables[0]]

        ## Performance test
        #for i in ['linearND', 'ndimage', 'nearest']:
        #    print i
        #    b = ReaderBlock(data, interpolation_method=i)
        #    preparation = b.interpolate(x, y)
        #    for num in range(5):
        #        t = datetime.now()
        #        values = b.interpolate_horizontal(slice2d)
        #        print num, np.mean(values[np.isfinite(values)]), datetime.now()-t

        block3d = np.ma.array([slice2d, 1.1*slice2d, 1.2*slice2d,
                               1.3*slice2d, 1.4*slice2d, 2*slice2d])

        b = ReaderBlock(data, interpolation_method='linearND')
        preparation = b.interpolate(x, y)
        #values = b.interpolate_horizontal_layers(slice2d)
        values = b.interpolate_horizontal_layers(block3d)
        print values.shape
        print values

        stop
        b = ReaderBlock(data, interpolation_method='linearND')
        preparation = b.interpolate(x, y)
        values = b.interpolate_horizontal(slice2d)

        b2 = ReaderBlock(data, interpolation_method='ndimage')
        preparation = b2.interpolate(x, y)
        values2 = b2.interpolate_horizontal(slice2d)

        plt.scatter(values, values2)
        plt.grid()
        plt.show()

        #self.assertEqual(covered.tolist(), [1, 2])
        #self.assertTrue(r.covers_time(r.start_time))
        #self.assertFalse(r.covers_time(r.start_time - r.time_step))
        #self.assertFalse(r.proj.is_latlong())


if __name__ == '__main__':
    unittest.main()
