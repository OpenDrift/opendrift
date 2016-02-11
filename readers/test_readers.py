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

    def test_vertical_interpolation_from_buffer(self):
        norkyst3d = reader_netCDF_CF_generic.Reader(script_folder +
            '/../test_data/14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        lon = np.array([4.73])
        lat = np.array([62.35])
        z = np.array([-33, 0])
        variables = ['x_sea_water_velocity', 'x_sea_water_velocity',
                     'sea_water_temperature']
        # Call get_variables_from_buffer which interpolates both in 
        # space (horizontally, vertically) and then in time
        data, profiles = norkyst3d.get_variables_from_buffer(
                variables, profiles='sea_water_temperature',
                profiles_depth = [-100, 0],
                time = norkyst3d.start_time + timedelta(seconds=900),
                lon=lon, lat=lat, z=z, block=True)
        # Check surface value
        self.assertEqual(data['sea_water_temperature'][1],
                         profiles['sea_water_temperature'][0])
        # Check interpolated temperature at 33 m depth
        self.assertAlmostEqual(data['sea_water_temperature'][0],
                               8.3167563152313)

    def atest_vertical_interpolation_sigma(self):
        norkyst3d = reader_ROMS_native.Reader(script_folder +
            '/../test_data/2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        # Default buffer may be changed in future, but this test case
        # requires that buffer is 25
        lon = np.array([12.46, 12.46])
        lat = np.array([68.21, 68.31])
        z = np.array([-33, 0])
        variables = ['x_sea_water_velocity', 'x_sea_water_velocity',
                     'sea_water_temperature']
        # Call get_variables_from_buffer which interpolates both in 
        # space (horizontally, vertically) and then in time
        data, profiles = norkyst3d.get_variables_from_buffer(
                variables, profiles='sea_water_temperature',
                profiles_depth = [-100, 0],
                time = norkyst3d.start_time + timedelta(seconds=900),
                lon=lon, lat=lat, z=z, block=True)
        #print data

        # Vertical interpolation of sigma to be implemented soon.
        # Currently z is simply neglected

        # Check surface value
        #print profiles
        #self.assertEqual(data['sea_water_temperature'][1],
        #                 profiles['sea_water_temperature'][0])
        ## Check interpolated temperature at 33 m depth
        #self.assertAlmostEqual(data['sea_water_temperature'][0],
        #                       8.3167563152313)


if __name__ == '__main__':
    unittest.main()
