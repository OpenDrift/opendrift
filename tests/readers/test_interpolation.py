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

import os
import unittest
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

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

    def test_dateline(self):

        # Make synthetic netCDF file with currents from 0 to 360 deg longitude
        fc = 'opendrift_test_current_0_360.nc'
        lon = np.arange(0, 360)
        lat = np.arange(-88, 89)
        start_time = datetime(2021, 1, 1)
        time = [start_time + i*timedelta(hours=24) for i in range(3)]
        t, xcurr, ycurr = np.meshgrid(time, np.zeros(lat.shape), np.zeros(lon.shape), indexing='ij')
        xcurr[:, :, 0:180] = 1  # eastward
        xcurr[:, :, 180:] = -1  # westward, i.e. current divergence at lon=0
        ds = xr.Dataset(
            {"xcurr": (("time", "lat", "lon"), xcurr, {'standard_name': 'x_sea_water_velocity'}),
             "ycurr": (("time", "lat", "lon"), ycurr, {'standard_name': 'y_sea_water_velocity'})},
            coords={"lon": lon, "lat": lat, "time": time})
        ds.to_netcdf(fc)

        # Make synthetic netCDF file with winds from -180 to 180 deg longitude
        fw = 'opendrift_test_winds_180_180.nc'
        lon = np.arange(-180, 180)
        t, xwind, ywind = np.meshgrid(time, np.zeros(lat.shape), np.zeros(lon.shape), indexing='ij')
        ywind[:, :, 0:180] = 1  # northward
        ywind[:, :, 180:] = -1  # southward, i.e. wind divergence at lon=180
        ds = xr.Dataset(
            {"xwind": (("time", "lat", "lon"), xwind, {'standard_name': 'x_wind'}),
             "ywind": (("time", "lat", "lon"), ywind, {'standard_name': 'y_wind'})},
            coords={"lon": lon, "lat": lat, "time": time})
        ds.to_netcdf(fw)

        # Make synthetic netCDF file with winds from 160 to 280 deg longitude (Pacific)
        fw2 = 'opendrift_test_winds_160_280.nc'
        lon = np.arange(160, 280)
        t, xwind, ywind = np.meshgrid(time, np.zeros(lat.shape), np.zeros(lon.shape), indexing='ij')
        ywind[:, :, 0:20] = -1  # southhward
        ywind[:, :, 20:] = 1  # northward, i.e. wind divergence at lon=180
        ds = xr.Dataset(
            {"xwind": (("time", "lat", "lon"), xwind, {'standard_name': 'x_wind'}),
             "ywind": (("time", "lat", "lon"), ywind, {'standard_name': 'y_wind'})},
            coords={"lon": lon, "lat": lat, "time": time})
        ds.to_netcdf(fw2)

        reader_current = reader_netCDF_CF_generic.Reader(fc)
        reader_wind = reader_netCDF_CF_generic.Reader(fw)
        reader_wind2 = reader_netCDF_CF_generic.Reader(fw2)

        lons = np.array([-175, 0, 175])
        lats = np.array([60, 60, 60])
        np.testing.assert_array_almost_equal(reader_wind.covers_positions(lons, lats)[0], [0, 1, 2], decimal=1)
        np.testing.assert_array_almost_equal(reader_wind2.covers_positions(lons, lats)[0], [0, 2], decimal=1)

        # Simulation across 0 meridian
        o = OceanDrift(loglevel=30)
        o.add_readers_from_list([fc, fw])
        o.seed_elements(lon=[-2, 2], lat=[60, 60], time=start_time, wind_drift_factor=.1)
        o.run(steps=2)
        # Check that current give divergence, and that
        # wind is northwards east of 0 and southwards to the east
        np.testing.assert_array_almost_equal(o.elements.lon, [-2.129,  2.129], decimal=3)
        np.testing.assert_array_almost_equal(o.elements.lat, [60.006, 59.994], decimal=3)

        # Simulation across dateline (180 E/W)
        o = OceanDrift(loglevel=30)
        o.add_readers_from_list([fc, fw])
        o.seed_elements(lon=[-175, 175], lat=[60, 60], time=start_time, wind_drift_factor=.1)
        o.run(steps=2)
        #o.plot(fast=True)
        # Check that current give convergence, and that
        # wind is northwards east of 180 and southwards to the west
        np.testing.assert_array_almost_equal(o.elements.lon, [-175.129,  175.129], decimal=3)
        np.testing.assert_array_almost_equal(o.elements.lat, [60.006, 59.994], decimal=3)

        # Same as above, but with wind reader from 160 to 280 deg
        o = OceanDrift(loglevel=30)
        o.add_readers_from_list([fc, fw2])
        o.seed_elements(lon=[-175, 175], lat=[60, 60], time=start_time, wind_drift_factor=.1)
        o.run(steps=2)
        #o.plot(fast=True)
        # Check that current give convergence, and that
        # wind is northwards east of 180 and southwards to the west
        np.testing.assert_array_almost_equal(o.elements.lon, [-175.129,  175.129], decimal=3)
        np.testing.assert_array_almost_equal(o.elements.lat, [60.006, 59.994], decimal=3)

        # Cleaning up
        os.remove(fw)
        os.remove(fw2)
        os.remove(fc)

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
        self.assertAlmostEqual(values[10], 1.6487979858538129)
        self.assertEqual(sum(values.mask), 15)

    def test_flipped(self):
        x = np.arange(10).astype(np.float32)
        y = np.arange(20).astype(np.float32)
        X, Y = np.meshgrid(x, y)
        d = X*5
        b = ReaderBlock({'x': x, 'y': y, 'v': d, 'time': datetime.now()})
        b_flipped = ReaderBlock({'x': np.flip(x), 'y': y, 'v': np.flip(d, axis=1), 'time': datetime.now()})
        x0 = np.array([0, 7, 7.3, 7.41, 9])  # Some random points
        y0 = np.array([5, 5, 8, 8.2, 5])
        bi = b.Interpolator2DClass(b.x, b.y, x0, y0)
        bi_flipped = b_flipped.Interpolator2DClass(b_flipped.x, b_flipped.y, x0, y0)

        np.testing.assert_array_almost_equal(bi(d), bi_flipped(np.flip(d, axis=1)))
        np.testing.assert_array_almost_equal(bi(d), 5*x0)

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
                                    x=x, y=y, z=z)

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
                                    x=x, y=y, z=z)

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
                                    x=x, y=y, z=z)

        # Introduce missing values
        data['x_sea_water_velocity'] = np.ma.masked_where(
            data['x_sea_water_velocity']>.08,
            data['x_sea_water_velocity'])

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
                                    x=x, y=y, z=z)

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
                                    x=x, y=y, z=0)
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
