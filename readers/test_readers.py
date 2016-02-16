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
import inspect
from datetime import timedelta, datetime

import numpy as np

import reader_netCDF_CF_generic
import reader_ROMS_native

script_folder = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe())))

class TestReaders(unittest.TestCase):
    """Tests for readers"""

    def test_reader_netcdf(self):
        """Check reader functionality."""

        reader1 = reader_netCDF_CF_generic.Reader(script_folder + 
            '/../test_data/16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        reader2 = reader_ROMS_native.Reader(script_folder +
            '/../test_data/2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        readers = [reader1, reader2]

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


    def test_vertical_profiles(self):

        norkyst3d = reader_netCDF_CF_generic.Reader(script_folder +
            '/../test_data/14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        lon = np.array([4.73])
        lat = np.array([62.35])
        variables = ['x_sea_water_velocity', 'x_sea_water_velocity',
                     'sea_water_temperature']
        x,y = norkyst3d.lonlat2xy(lon, lat)
        data = norkyst3d.get_variables(variables,
                                       time=norkyst3d.start_time,
                                       x=x, y=y, z=[0, -100], block=True)
        self.assertEqual(data['z'][4], -25)
        self.assertEqual(data['z'][4], -25)
        self.assertAlmostEqual(data['sea_water_temperature'][:,0,0][7],
                         9.220000267028809)

    def test_vertical_interpolation(self):
        norkyst3d = reader_netCDF_CF_generic.Reader(script_folder +
            '/../test_data/14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        lon = np.array([4.73, 4.75])
        lat = np.array([62.35, 62.30])
        z = np.array([0, -33])
        variables = ['x_sea_water_velocity', 'x_sea_water_velocity',
                     'sea_water_temperature']
        # Call get_variables_interpolated which interpolates both in 
        # space (horizontally, vertically) and then in time
        data, profiles = norkyst3d.get_variables_interpolated(
                variables, profiles=['sea_water_temperature'],
                profiles_depth = [-100, 0],
                time = norkyst3d.start_time + timedelta(seconds=900),
                lon=lon, lat=lat, z=z, block=True)
        # Check surface value
        self.assertEqual(data['sea_water_temperature'][0],
                         profiles['sea_water_temperature'][0,0])
        # Check interpolated temperature at 33 m depth
        self.assertAlmostEqual(data['sea_water_temperature'][1],
                               8.2648999309539786)
        #import matplotlib.pyplot as plt
        #plt.plot(profiles['sea_water_temperature'][:,0])
        #plt.plot(profiles['sea_water_temperature'][:,1], 'r')
        #plt.show()

    def test_vertical_interpolation_sigma(self):
        nordic3d = reader_ROMS_native.Reader(script_folder +
            '/../test_data/2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        lon = np.array([12.46, 12.46, 12.46])
        lat = np.array([68.21, 69.31, 69.31])
        z = np.array([-33, 0, -2500])
        x, y = nordic3d.lonlat2xy(lon, lat)
        variables = ['x_sea_water_velocity', 'x_sea_water_velocity',
                     'sea_water_temperature']
        # Call get_variables_interpolated which interpolates both in 
        data = nordic3d.get_variables(variables,
                time = nordic3d.start_time + timedelta(seconds=900),
                x=x, y=y, z=z, block=True)
        #print data
        self.assertAlmostEqual(data['sea_water_temperature'][0,60, 60],
                               3.4470012188)
        self.assertAlmostEqual(data['sea_water_temperature'][-1,60, 60],
                               -0.297809013466)
        
if __name__ == '__main__':
    unittest.main()
