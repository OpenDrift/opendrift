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
from opendrift.models.leeway import Leeway
from opendrift.models.openoil3D import OpenOil3D
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_constant
from opendrift.readers import reader_lazy
from opendrift.readers import reader_from_url
from opendrift.models.pelagicegg import PelagicEggDrift


o = OceanDrift(loglevel=20)

reader_list = [
    'www.nonexistingurl.com',
    o.test_data_folder() +
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc',
    '/nonexistingdisk/nonexistingfile.ext',
    o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF.nc_20160202_subset']


class TestReaders(unittest.TestCase):
    """Tests for readers"""

    def test_adding_readers(self):
        o = OceanDrift()
        basemap = reader_basemap_landmask.Reader(
            llcrnrlon=-1.5, llcrnrlat=59,
            urcrnrlon=7, urcrnrlat=64, resolution='c')
        r = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        o.add_reader([r, basemap])
        self.assertEqual(o.priority_list['land_binary_mask'],
                         ['roms native', 'basemap_landmask'])
        self.assertEqual(o.priority_list['x_sea_water_velocity'],
                         ['roms native'])
        # Switch order
        o = OceanDrift()
        o.add_reader([basemap, r])
        self.assertEqual(o.priority_list['land_binary_mask'],
                         ['basemap_landmask', 'roms native'])
        self.assertEqual(o.priority_list['x_sea_water_velocity'],
                         ['roms native'])

        # Test add_readers_from_list
        o = OceanDrift()
        o.add_readers_from_list(reader_list, lazy=False)
        self.assertEqual(o.priority_list['x_sea_water_velocity'],
                         ['roms native'])
        self.assertEqual(o.priority_list['x_wind'],
                         [o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF.nc_20160202_subset'])

    def test_reader_from_url(self):
        readers = reader_from_url(reader_list)
        self.assertIsNone(readers[0])
        self.assertTrue(isinstance(readers[1],
                                   reader_ROMS_native.Reader))
        self.assertIsNone(readers[2])
        self.assertTrue(isinstance(readers[3],
                                   reader_netCDF_CF_generic.Reader))

    def test_lazy_reader(self):
        o = OceanDrift(loglevel=20)
        lr = reader_lazy.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        self.assertFalse(lr.initialised)
        self.assertEqual(len(lr.covers_positions([15], [69])), 1)
        self.assertEqual(len(lr.covers_positions([0], [0])), 0)
        self.assertTrue(lr.initialised)

        # Make a corresponding, unlazy reader
        rr = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        self.assertEqual(len(rr.covers_positions([15], [69])), 1)
        self.assertEqual(len(rr.covers_positions([0], [0])), 0)

        # Check that both readers provide the same attributes
        for att in rr.__dict__:
            self.assertEqual(type(lr.__getattr__(att)),
                             type(getattr(rr, att)))
            if type(getattr(rr, att)) in [float, int, dict, str, list,
                                          datetime, timedelta, bool,
                                          np.float64]:
                self.assertEqual(lr.__getattr__(att),
                                 getattr(rr, att))
            elif type(getattr(rr, att)) in [np.ndarray]:
                self.assertIsNone(np.testing.assert_array_equal(
                                  lr.__getattr__(att),
                                  getattr(rr, att)))
            else:
                print('Skipping: ' + att + ' ' +
                      str(type(getattr(rr, att))))

    def test_lazy_reader_oildrift(self):
        o = OpenOil3D(loglevel=20)
        reader_constant_wind = \
            reader_constant.Reader({'x_wind':5, 'y_wind': 6})
        o.add_reader(reader_constant_wind)

        o.add_readers_from_list(reader_list, lazy=True)

        self.assertEqual(len(o._lazy_readers()), 4)
        o.set_config('general:basemap_resolution', 'c')
        o.seed_elements(lon=14, lat=67.85,
                        time=datetime(2016, 2, 2, 12))
        o.run(steps=5)
        self.assertEqual(len(o._lazy_readers()), 2)
        self.assertEqual(len(o.discarded_readers), 1)

    #def test_oildrift_backwards(self):
    #    o = OpenOil3D(loglevel=20)
    #    reader_constant_wind = \
    #        reader_constant.Reader({'x_wind':5, 'y_wind': 6})
    #    o.add_reader(reader_constant_wind)

    #    o.add_readers_from_list(reader_list, lazy=True)

    #    self.assertEqual(len(o._lazy_readers()), 4)
    #    o.set_config('general:basemap_resolution', 'c')
    #    o.seed_elements(lon=14, lat=67.85,
    #                    time=datetime(2016, 2, 2, 12))
    #    o.set_config()
    #    o.run(steps=5)
    #    self.assertEqual(len(o._lazy_readers()), 2)
    #    self.assertEqual(len(o.discarded_readers), 1)

    #def test_lazy_reader_oildrift_real(self):
    #    o = OpenOil3D(loglevel=0)
    #    o.add_readers_from_file(o.test_data_folder() +
    #        '../../opendrift/scripts/data_sources.txt')

    #    o.set_config('general:basemap_resolution', 'c')
    #    o.seed_elements(lon=4, lat=60.0,
    #                    time=datetime(2018, 7, 2, 12))
    #    o.run(steps=5)
    #    print o

    def test_lazy_reader_leeway_compare(self):

        o1 = Leeway(loglevel=0)
        #o1.fallback_values['land_binary_mask'] = 0
        o1.required_variables = [r for r in o1.required_variables
                                 if r != 'land_binary_mask']
        o1.add_readers_from_list(reader_list, lazy=False)
        time = o1.readers['roms native'].start_time
        o1.seed_elements(lat=67.85, lon=14, time=time)
        o1.run(steps=5)

        o2 = Leeway(loglevel=20)
        #o2.fallback_values['land_binary_mask'] = 0
        o2.required_variables = [r for r in o1.required_variables
                                 if r != 'land_binary_mask']
        o2.add_readers_from_list(reader_list, lazy=True)
        o2.seed_elements(lat=67.85, lon=14, time=time)
        o2.run(steps=5)

        # Some differences in wind and current components
        # due to different coordinate system
        for var in o1.history.dtype.names:
            if var in ['x_wind', 'y_wind', 'x_sea_water_velocity',
                       'y_sea_water_velocity']:
                tolerance = 1
            else:
                tolerance = 5
            self.assertIsNone(np.testing.assert_array_almost_equal(
                o1.history[var], o2.history[var], tolerance))

    def test_constant_and_lazy_reader_leeway(self):
        cw = reader_constant.Reader({'x_wind':5, 'y_wind': 6})
        cc = reader_constant.Reader({'x_sea_water_velocity':0,
                                     'y_sea_water_velocity': .2})

        o = Leeway(loglevel=20)
        o.set_config('general:basemap_resolution', 'c')
        o.add_reader([cw, cc])
        o.add_readers_from_list(reader_list)
        o.fallback_values['x_sea_water_velocity'] = 0.0
        o.fallback_values['y_sea_water_velocity'] = 0.1
        time = datetime(2016,2,2,12)
        o.seed_elements(lat=67.85, lon=14, time=time)
        o.run(steps=2)
        self.assertAlmostEqual(o.elements.lat[0], 67.8791, 3)

    def test_automatic_basemap(self):
        self.assertRaises(ValueError, o.run)
        o.seed_elements(lon=4, lat=60, time=datetime(2016,9,1))
        o.set_config('general:basemap_resolution', 'c')  # To make test fast
        o.run(steps=2)

    def test_reader_coverage(self):
        r = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        # Element outside reader domain
        self.assertEqual(len(r.covers_positions(5, 80)), 0)
        x, y = r.lonlat2xy(5, 80)
        self.assertRaises(ValueError, r.check_arguments,
                          'y_sea_water_velocity', r.start_time, x, y, 0)
        # Element inside reader domain
        self.assertEqual(len(r.covers_positions(5, 60)), 1)
        x, y = r.lonlat2xy(5, 60)
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
            print(r)
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
                    print('#'*60)
                    print('#'*60)
                    print('WARNING: A point on the boundary is considered ' \
                          'outside after conversion x,y -> lon,lat -> x,y. ' \
                          'This is different from "standard", but is due to ' \
                          'rounding differences and not considered to be an ' \
                          'error. Numpy version is %s' % (np.__version__))
                    print('#'*60)
                    print('#'*60)
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
                               3.447, 2)
                               #3.59, 2)
        self.assertAlmostEqual(data['sea_water_temperature'][-1,60, 60],
                               -0.783, 2)
                               #-0.803, 2)

    def test_get_environment(self):
        o = PelagicEggDrift(loglevel=30)
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
        self.assertIsNone(np.testing.assert_array_almost_equal(
            [0], reader_nordic.covers_positions(testlon, testlat, testz)))
        self.assertIsNone(np.testing.assert_array_almost_equal(
            [0, 1, 2],
            reader_arctic.covers_positions(testlon, testlat, testz)))
        o.seed_elements(testlon, testlat, testz, time=reader_nordic.start_time)
        o.fallback_values['land_binary_mask'] = 0
        env, env_profiles, missing = \
            o.get_environment(o.required_variables,
                              reader_nordic.start_time,
                              testlon, testlat, testz,
                              o.required_profiles)
        self.assertAlmostEqual(env['sea_water_temperature'][0], 4.267, 2)
        self.assertAlmostEqual(env['sea_water_temperature'][1], 0.468122, 3)
        self.assertAlmostEqual(env['sea_water_temperature'][4], 10.0)
        self.assertIsNone(np.testing.assert_array_almost_equal(
            missing, [False,False,False,False,False]))
        self.assertAlmostEqual(env_profiles['sea_water_temperature'][0,0],
                               4.267, 2)
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
        self.assertEqual(set(env_profiles2.keys()),
            set(['z', 'sea_water_temperature']))
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

        self.assertIsNone(np.testing.assert_array_almost_equal(
            env['x_sea_water_velocity'],
            env2['x_sea_water_velocity']))
        self.assertIsNone(np.testing.assert_array_almost_equal(
            env_profiles2['sea_water_temperature'].ravel(),
            env_profiles4['sea_water_temperature'].ravel()))


    def test_constant_reader(self):
        o = OpenOil3D(loglevel=0)
        o.set_config('general:basemap_resolution', 'c')
        cw = reader_constant.Reader({'x_wind':5, 'y_wind': 6})
        cc = reader_constant.Reader({'x_sea_water_velocity':0, 'y_sea_water_velocity': .2})
        cs = reader_constant.Reader({'sea_water_temperature': 278})
        r = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        o.add_reader([cw, cc, r])
        # TODO: should check why adding constant reader with 
        #   sea_water_temperature gives Deprecated warning
        #o.add_reader([cw, cc, cs, r])
        o.seed_elements(lon=4, lat=60, time=r.start_time, number=5)
        o.run(steps=3)

if __name__ == '__main__':
    unittest.main()
