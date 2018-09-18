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
from datetime import datetime

import numpy as np

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers.interpolation import \
        expand_numpy_array, \
        ReaderBlock, LinearND2DInterpolator, \
        NDImage2DInterpolator, Nearest2DInterpolator, \
        Nearest1DInterpolator, Linear1DInterpolator

o = OceanDrift()


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
        z = np.linspace(data_dict['z'].min(), data_dict['z'].max(), 100)

        return data_dict, x, y, z

    def test_covers_positions(self):
        
        data_dict, x, y, z = self.get_synthetic_data_dict()
        # Make block from dictionary, and apply tests
        b = ReaderBlock(data_dict)

        xn = np.linspace(-70, 470, 100)
        yn = np.linspace(10, 340, 100)
        self.assertTrue(b.covers_positions(xn, yn))

        xn = np.linspace(500, 600, 100)
        yn = np.linspace(400, 500, 100)
        self.assertFalse(b.covers_positions(xn, yn))

        xn = np.linspace(400, 500, 100)
        yn = np.linspace(0, 30, 100)
        self.assertFalse(b.covers_positions(xn, yn))

    def test_interpolation_horizontal(self):

        data_dict, x, y, z = self.get_synthetic_data_dict()
        # Make block from dictionary, and apply tests
        b = ReaderBlock(data_dict, interpolation_horizontal='ndimage')
        self.assertEqual(b.data_dict['var2d'].shape,
                         (len(b.y), len(b.x)))
        self.assertEqual(b.data_dict['var3d'].shape,
                         (len(b.z), len(b.y), len(b.x)))
        # Make some element positions
        interpolator2d = b.Interpolator2DClass(b.x, b.y, x, y)
        values = interpolator2d(data_dict['var2d'])
        # Checking output is as expected
        self.assertEqual(values[10], 1.6487979858538129)
        self.assertEqual(sum(values.mask), 15)

    def test_interpolation_ensemble(self):
        
        data_dict, x, y, z = self.get_synthetic_data_dict()
        x = x[0:15]
        y = y[0:15]
        z = z[0:15]
        data_dict['var2d'] = np.ones(data_dict['var2d'].shape)
        data_dict['var3d'] = np.ones(data_dict['var3d'].shape)
        data_dict['var2de'] = [data_dict['var2d']*1,
                               data_dict['var2d']*2,
                               data_dict['var2d']*3]
        data_dict['var3de'] = [data_dict['var3d']*31,
                               data_dict['var3d']*32,
                               data_dict['var3d']*33]

        b = ReaderBlock(data_dict)
        interp = b.interpolate(x, y, z)[0]  # 1 is profiles
        v2 = interp['var2d']
        v2e = interp['var2de']
        v3 = interp['var3d']
        v3e = interp['var3de']
        
        self.assertEqual(v2[0], 1)
        self.assertEqual(v2e[0], 1)
        self.assertEqual(v2e[1], 2)
        self.assertEqual(v2e[3], 1)
        self.assertEqual(v3[0], 1)
        self.assertEqual(v3e[0], 31)
        self.assertEqual(v3e[1], 32)
        self.assertAlmostEqual(v3e[3], 31)

    def test_interpolation_vertical(self):

        # 3 elements, 4 depths
        zgrid = np.array([0, 1, 3, 10])
        z = np.array([.5, 3, 9])
        data = np.array([[0, 0, 0],
                         [1, 1, 1],
                         [2, 2, 2],
                         [3, 3, 3]])
        interpolator = Nearest1DInterpolator(zgrid, z)
        self.assertTrue(np.allclose(interpolator(data), [0, 2, 3]))
        interpolator = Linear1DInterpolator(zgrid, z)
        self.assertTrue(np.allclose(interpolator(data),
                                    [0.5, 2, 2.85714286]))

        # And with exatrapolation (~to surface and bottom)
        zgrid = np.array([1, 3, 5, 10])
        z = np.array([.5, 6, 12])
        interpolator = Nearest1DInterpolator(zgrid, z)
        self.assertTrue(np.allclose(interpolator(data), [0, 2, 3]))
        interpolator = Linear1DInterpolator(zgrid, z)
        self.assertTrue(np.allclose(interpolator(data),
                                    [0.0, 2.2, 3]))

    def test_compare_interpolators(self):

        data_dict, x, y, z = self.get_synthetic_data_dict()
        arr = data_dict['var2d']
        # Make block from dictionary, and apply tests
        linearData = LinearND2DInterpolator(data_dict['x'], data_dict['y'],
                                            x, y)(data_dict['var2d'])
        nearestData = Nearest2DInterpolator(data_dict['x'], data_dict['y'],
                                            x, y)(data_dict['var2d'])
        ndimageData = NDImage2DInterpolator(data_dict['x'], data_dict['y'],
                                            x, y)(data_dict['var2d'])

        # Check that all interpolator give nearly equal values
        # for a given position
        i = 10
        self.assertAlmostEqual(linearData[i], nearestData[i], places=2)
        self.assertAlmostEqual(linearData[i], ndimageData[i], places=2)

    def test_interpolation_3dArrays(self):
        """Test interpolation."""
        reader = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
            '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')

        # 100000 points within 50x50 pixels over sea (corner of domain)
        num_points = 1000
        np.random.seed(0)  # To get the same random numbers each time
        x = np.random.uniform(reader.xmin, reader.xmin+800*50, num_points)
        y = np.random.uniform(reader.ymax-800*50, reader.ymax, num_points)
        z = np.random.uniform(-200, 0, num_points)
        variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                     'sea_water_temperature']
        # Read a block of data covering the points
        data = reader.get_variables(variables, time=reader.start_time,
                                    x=x, y=y, z=z, block=True)

        b = ReaderBlock(data, interpolation_horizontal='nearest')
        env, prof = b.interpolate(x, y, z, variables,
                                  profiles=['sea_water_temperature'],
                                  profiles_depth=[-30, 0])
        self.assertAlmostEqual(env['x_sea_water_velocity'][100],
                               0.075019, 3)
        self.assertAlmostEqual(prof['sea_water_temperature'][0,11],
                               7.549999, 3)
        self.assertAlmostEqual(prof['sea_water_temperature'][-1,11],
                               8.389999, 3)
        self.assertEqual(prof['z'][-1], b.z[-1])

    def test_zNone(self):
        d = {}
        d['x'] = np.arange(5)
        d['y'] = np.arange(7)
        d['z'] = 0
        d['time'] = None
        d['var'] = np.random.rand(5,6)
        rb = ReaderBlock(d)
        z = None
        i,p = rb.interpolate(np.array([1, 2]), np.array([2, 3]), z, 'var')
        self.assertTrue(i['var'][0] > 0)

    def test_repeated(self):
        """Check that block can be used for interpolation to several sets of positions"""
        reader = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
            '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')

        # 100 points within 50x50 pixels over sea (corner of domain)
        num_points = 100
        np.random.seed(0)  # To get the same random numbers each time
        x = np.random.uniform(reader.xmin, reader.xmin+800*50, num_points)
        y = np.random.uniform(reader.ymax-800*50, reader.ymax, num_points)
        z = np.random.uniform(-200, 0, num_points)
        variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                     'sea_water_temperature']
        # Read a block of data covering the points
        data = reader.get_variables(variables, time=reader.start_time,
                                    x=x, y=y, z=z, block=True)

        b = ReaderBlock(data, interpolation_horizontal='nearest')

        env, prof = b.interpolate(x, y, z, 'sea_water_temperature')
        x2 = x[20:30]
        y2 = y[20:30]
        z2 = z[20:30]
        env2, prof2 = b.interpolate(x2, y2, z2, 'sea_water_temperature')
        env3, prof3 = b.interpolate(x, y, z, 'sea_water_temperature')
        self.assertEqual(env['sea_water_temperature'][0],
                         env3['sea_water_temperature'][0])
        self.assertEqual(env['sea_water_temperature'][20],
                         env2['sea_water_temperature'][0])


    def test_interpolation_missing(self):
        """Test interpolation."""
        reader = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        num_points = 50
        np.random.seed(0)  # To get the same random numbers each time
        lons = np.random.uniform(10, 11, num_points)
        lats = np.random.uniform(66, 67.0, num_points)
        z = np.random.uniform(-200, 0, num_points)
        x, y = reader.lonlat2xy(lons, lats)

        variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                     'sea_water_temperature']
        # Read a block of data covering the points
        data = reader.get_variables(variables, time=reader.start_time,
                                    x=x, y=y, z=z, block=True)

        # Introduce missing values
        data['x_sea_water_velocity'].mask[[data['x_sea_water_velocity']>.08]] = True

        b = ReaderBlock(data, interpolation_horizontal='linearND')

        env, prof = b.interpolate(x, y, z, variables,
                                  profiles=['x_sea_water_velocity'],
                                  profiles_depth=[-30, 0])
        self.assertAlmostEqual(env['x_sea_water_velocity'][10],
                               0.074, 2)
        self.assertAlmostEqual(prof['x_sea_water_velocity'][5,48],
                               -0.090, 2)

    def test_linearNDFast(self):
        """Test interpolation."""
        reader = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        reader.buffer=3
        num_points = 50
        np.random.seed(0)  # To get the same random numbers each time
        lons = np.random.uniform(14, 15, num_points)
        lats = np.random.uniform(68, 68.4, num_points)
        z = np.random.uniform(-20, 0, num_points)
        x, y = reader.lonlat2xy(lons, lats)
        #z=None

        variables = ['x_sea_water_velocity']
        # Read a block of data covering the points
        data = reader.get_variables(variables, time=reader.start_time,
                                    x=x, y=y, z=z, block=True)

        b = ReaderBlock(data.copy(),
                        interpolation_horizontal='linearNDFast')
        x2 = np.random.uniform(x.min(), x.max(), num_points)
        y2 = np.random.uniform(y.min(), y.max(), num_points)
        z2 = np.random.uniform(-20, 0, num_points)
        self.assertTrue(b.covers_positions(x, y, z))
        self.assertTrue(b.covers_positions(x2, y2, z2))
        # Check that there are holes in the arrays of the ReaderBlock
        self.assertEqual(
            np.sum(~np.isfinite(b.data_dict['x_sea_water_velocity'])), 1001)
        # Check that LinearNDFast interpolation gives a real value
        env, prof = b.interpolate(x2, y2, z2,
                                  variables,
                                  profiles=variables,
                                  profiles_depth=[-30, 0])
        self.assertEqual(
            np.sum(~np.isfinite(env['x_sea_water_velocity'])), 0)
        # Check that the arrays of the ReaderBlock have been filled in
        self.assertEqual(
            np.sum(~np.isfinite(b.data_dict['x_sea_water_velocity'])), 0)

        # Check that nearest interpolation contains some NaN values
        b2 = ReaderBlock(data.copy(), interpolation_horizontal='nearest')
        env, prof = b2.interpolate(x2, y2, z2,
                                  variables,
                                  profiles=variables,
                                  profiles_depth=[-30, 0])
        self.assertEqual(
            np.sum(~np.isfinite(env['x_sea_water_velocity'])), 31)

    def test_expand_array(self):
        reader = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        reader.buffer=1
        num_points = 50
        np.random.seed(0)  # To get the same random numbers each time
        lons = np.random.uniform(14, 15, num_points)
        lats = np.random.uniform(68, 68.4, num_points)
        x, y = reader.lonlat2xy(lons, lats)
        variables = ['x_sea_water_velocity']
        # Read a block of data covering the points
        data = reader.get_variables(variables, time=reader.start_time,
                                    x=x, y=y, z=0, block=True)
        data = np.ma.filled(data['x_sea_water_velocity'],
                            fill_value=np.nan)
        self.assertTrue(np.isnan(data.max()))
        self.assertEqual(sum(~np.isfinite(data.ravel())), 80)
        expand_numpy_array(data)
        self.assertEqual(sum(~np.isfinite(data.ravel())), 40)
        expand_numpy_array(data)
        self.assertEqual(sum(~np.isfinite(data.ravel())), 9)
        expand_numpy_array(data)
        self.assertEqual(sum(~np.isfinite(data.ravel())), 0)
        self.assertFalse(np.isnan(data.max()))


if __name__ == '__main__':
    unittest.main()
