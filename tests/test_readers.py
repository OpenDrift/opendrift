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

from opendrift.models.oceandrift import OceanDrift
from opendrift.models.openoil3D import OpenOil3D
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_constant
from opendrift.models.pelagicegg import PelagicEggDrift


o = OceanDrift()
basemap = reader_basemap_landmask.Reader(
            llcrnrlon=-1.5, llcrnrlat=59,
            urcrnrlon=7, urcrnrlat=64, resolution='c')

class TestReaders(unittest.TestCase):
    """Tests for readers"""

    def test_adding_readers(self):
        o = OceanDrift()
        r = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        o.add_reader([r, basemap])
        #self.assertEqual(o.priority_list['land_binary_mask'],
        #                 ['roms native', 'basemap_landmask'])
        self.assertEqual(o.priority_list['x_sea_water_velocity'],
                         ['roms native'])
        # Switch order
        o = OceanDrift()
        o.add_reader([basemap, r])
        #self.assertEqual(o.priority_list['land_binary_mask'],
        #                 ['basemap_landmask', 'roms native'])
        self.assertEqual(o.priority_list['x_sea_water_velocity'],
                         ['roms native'])

    def test_automatic_basemap(self):
        self.assertRaises(ValueError, o.run)
        o.seed_elements(lon=4, lat=60, time=datetime(2016,9,1))
        o.set_config('general:basemap_resolution', 'c')  # To make test fast
        o.run(steps=2)

    def test_reader_coverage(self):
        r = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        # Element outside reader domain
        self.assertRaises(ValueError, r.check_coverage, r.start_time, 5, 80)
        x, y = r.lonlat2xy(5, 80)
        self.assertRaises(ValueError, r.check_arguments,
                          'y_sea_water_velocity', r.start_time, x, y, 0)
        # Element inside reader domain
        x, y, ind = r.check_coverage(r.start_time, 5, 60)  # inside
        self.assertEqual(ind, 0)
        var, time, x2, y2, z2, outside = \
            r.check_arguments('y_sea_water_velocity', r.start_time, x, y, 0)
        self.assertEqual(var, ['y_sea_water_velocity'])
        self.assertEqual(time, r.start_time)
        self.assertEqual(x, x2)
        self.assertEqual(y, y2)
        self.assertEqual(0, z2)
        self.assertEqual(len(outside), 0)

    def test_outside_reader_time_coverage(self):
        o = PelagicEggDrift()
        reader = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        o.add_reader(reader)
        o.fallback_values['x_sea_water_velocity'] = 1
        o.fallback_values['land_binary_mask'] = 0
        o.set_config('processes:turbulentmixing', False)
        o.seed_elements(lon=4.8, lat=60, number=1, time=reader.end_time)
        o.run(steps=2)
        # Check that fallback value is used when outside time coverage
        self.assertEqual(o.history['x_sea_water_velocity'][0][-1], 1.0)

    def test_reader_netcdf(self):
        """Check reader functionality."""

        reader1 = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        reader2 = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
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
            if len(covered) != 1:
                self.assertEqual(covered.tolist(), [1, 2])
            else:
                if covered == [2]:
                    print '#'*60
                    print '#'*60
                    print 'WARNING: A point on the boundary is considered ' \
                          'outside after conversion x,y -> lon,lat -> x,y. ' \
                          'This is different from "standard", but is due to ' \
                          'rounding differences and not considered to be an ' \
                          'error. Numpy version is %s' % (np.__version__)
                    print '#'*60
                    print '#'*60
                else:
                    self.assertTrue(False)  # Should never happen!

            self.assertTrue(r.covers_time(r.start_time))
            self.assertFalse(r.covers_time(r.start_time - r.time_step))
            self.assertFalse(r.proj.is_latlong())


    def test_vertical_profiles(self):

        norkyst3d = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
            '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
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
        norkyst3d = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
            '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
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
                               8.36, 2)
        #import matplotlib.pyplot as plt
        #plt.plot(profiles['sea_water_temperature'][:,0])
        #plt.plot(profiles['sea_water_temperature'][:,1], 'r')
        #plt.show()

    def test_vertical_interpolation_sigma(self):
        nordic3d = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        lon = np.array([12.46, 12.46, 12.46])
        lat = np.array([68.21, 69.31, 69.31])
        z = np.array([-33, 0, -2500])
        x, y = nordic3d.lonlat2xy(lon, lat)
        variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                     'sea_water_temperature']
        # Call get_variables_interpolated which interpolates both in 
        data = nordic3d.get_variables(variables,
                time = nordic3d.start_time + timedelta(seconds=900),
                x=x, y=y, z=z, block=True)
        self.assertAlmostEqual(data['sea_water_temperature'][0,60, 60],
                               3.59, 2)
        self.assertAlmostEqual(data['sea_water_temperature'][-1,60, 60],
                               -0.803, 2)

    def test_get_environment(self):
        o = PelagicEggDrift(loglevel=0)
        reader_nordic = reader_ROMS_native.Reader(o.test_data_folder() + '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc', name='Nordic')
        reader_arctic = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '2Feb2016_Nordic_sigma_3d/Arctic20_1to5Feb_2016.nc', name='Arctic')
        ######################################################
        # Vertical interpolation is another issue to be fixed:
        reader_nordic.zlevels = reader_arctic.z
        ######################################################
        o.add_reader([reader_nordic, reader_arctic])
        # One point covered only by Nordic, two points coverd
        # by both readers, and two points covered by none of the readers
        testlon = np.array((14.0, 20.0, 20.1, 4, 5))
        testlat = np.array((70.1, 76.0, 76.1, 60, 60))
        testz = np.random.uniform(0, 0, len(testlon))
        self.assertItemsEqual([0], reader_nordic.covers_positions(
                                    testlon, testlat, testz))
        self.assertItemsEqual([0, 1, 2], reader_arctic.covers_positions(
                                    testlon, testlat, testz))
        o.seed_elements(testlon, testlat, testz, time=reader_nordic.start_time)
        o.fallback_values['land_binary_mask'] = 0
        env, env_profiles, missing = \
            o.get_environment(o.required_variables,
                              reader_nordic.start_time,
                              testlon, testlat, testz,
                              o.required_profiles)
        self.assertAlmostEqual(env['sea_water_temperature'][0], 4.318, 2)
        self.assertAlmostEqual(env['sea_water_temperature'][1], 0.468122, 3)
        self.assertAlmostEqual(env['sea_water_temperature'][4], 10.0)
        self.assertItemsEqual(missing, [False,False,False,False,False])
        self.assertAlmostEqual(env_profiles['sea_water_temperature'][0,0],
                               4.318, 2)
        self.assertAlmostEqual(env_profiles['sea_water_temperature'][0,4], 10)
        self.assertAlmostEqual(env_profiles['sea_water_temperature'][8,2], 10)
        self.assertAlmostEqual(env_profiles['sea_water_temperature'][7,2],
                               2.3049809, 3)
        # Get separate data
        env2, env_profiles2, missing2 = \
            o.get_environment(['x_sea_water_velocity', 'y_sea_water_velocity',
                               'sea_water_temperature'],
                              reader_nordic.start_time,
                              testlon, testlat, testz,
                              ['sea_water_temperature'])
        self.assertTrue(env_profiles2 is not None)
        self.assertEqual(env_profiles2.keys(), ['z', 'sea_water_temperature'])
        # Get separate data, without profile
        env3, env_profiles3, missing3 = \
            o.get_environment(['x_sea_water_velocity', 'y_sea_water_velocity',
                               'sea_water_temperature'],
                              reader_nordic.start_time,
                              testlon, testlat, testz,
                              profiles=None)
        self.assertTrue(env_profiles3 is None)
        # Get separate data
        env4, env_profiles4, missing4 = \
            o.get_environment(['x_sea_water_velocity', 'y_sea_water_velocity',
                               'sea_water_temperature'],
                              reader_nordic.start_time,
                              testlon, testlat, testz,
                              ['sea_water_temperature'])

        self.assertItemsEqual(env['x_sea_water_velocity'],
                              env2['x_sea_water_velocity'])
        #print env_profiles['sea_water_temperature'], '1'*50
        #print env_profiles2['sea_water_temperature'], '2'*50
        #print env_profiles4['sea_water_temperature'], '4'*50
        # Test below should also pass, To be fixed
        #self.assertItemsEqual(env_profiles['sea_water_temperature'].ravel(),
        #                      env_profiles2['sea_water_temperature'].ravel())
        self.assertItemsEqual(env_profiles2['sea_water_temperature'].ravel(),
                              env_profiles4['sea_water_temperature'].ravel())


    def test_constant_reader(self):
        o = OpenOil3D(loglevel=0)
        o.set_config('general:basemap_resolution', 'c')
        cw = reader_constant.Reader({'x_wind':5, 'y_wind': 6})
        cc = reader_constant.Reader({'x_sea_water_velocity':0, 'y_sea_water_velocity': .2})
        cs = reader_constant.Reader({'sea_water_temperature': 278})
        r = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        o.add_reader([cw, cc, cs, r])
        o.seed_elements(lon=4, lat=60, time=r.start_time, number=5)
        o.run(steps=3)

if __name__ == '__main__':
    unittest.main()
