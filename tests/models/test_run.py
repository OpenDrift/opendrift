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

import unittest
import pytest
from datetime import datetime, timedelta
import os
import inspect

import numpy as np

from opendrift.readers import reader_ArtificialOceanEddy
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_oscillating
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.basemodel import Mode, WrongMode
from opendrift.models.openoil import OpenOil
from opendrift.models.leeway import Leeway
from opendrift.models.pelagicegg import PelagicEggDrift
from opendrift.models.plastdrift import PlastDrift


def gdal_error_handler(err_class, err_num, err_msg):
    errtype = {
            gdal.CE_None:'None',
            gdal.CE_Debug:'Debug',
            gdal.CE_Warning:'Warning',
            gdal.CE_Failure:'Failure',
            gdal.CE_Fatal:'Fatal'
    }
    err_msg = err_msg.replace('\n',' ')
    err_class = errtype.get(err_class, 'None')
    print('Error Number: %s' % (err_num))
    print('Error Type: %s' % (err_class))
    print('Error Message: %s' % (err_msg))

try:
    from osgeo import ogr, osr, gdal
    version_num = int(gdal.VersionInfo('VERSION_NUM'))
    if version_num >= 2000000:
        has_ogr = True
        gdal.PushErrorHandler(gdal_error_handler)
    else:
        print('GDAL version >= 2.0 is required:')
        has_ogr = False
except Exception as e:
    print('GDAL is not available:')
    print(e)
    has_ogr = False

class TestRun(unittest.TestCase):
    """Tests for (non-scalar) LagrangianArray"""

    def make_OceanDrift_object(self):
        self.o = OceanDrift(loglevel=30)
        self.fake_eddy = reader_ArtificialOceanEddy.Reader(2, 62)
        self.reader_landmask = reader_global_landmask.Reader()
        self.o.add_reader([self.fake_eddy, self.reader_landmask])

    def test_seed(self):
        """Test seeding"""
        o = OceanDrift(loglevel=20)
        number = 3
        lonvec = np.linspace(2, 5, number)
        latvec = np.linspace(60, 61, number)
        o.seed_elements(lonvec, latvec, number=number,
                        time=datetime(2015, 1, 1, 12, 5, 17))
                        #time=[datetime(2015, 1, 1), datetime(2015, 1, 3)])

        # Check that 6 elements are scheduled, but none seeded
        self.assertEqual(o.num_elements_scheduled(), number)
        self.assertEqual(o.num_elements_active(), 0)
        self.assertEqual(o.num_elements_activated(), 0)
        self.assertEqual(o.num_elements_deactivated(), 0)
        self.assertEqual(o.num_elements_total(), number)

    def test_modes(self):
        '''Test that methods are not allowed if wrong mode'''
        o = OceanDrift(loglevel=50)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
                    '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o.set_config('environment:constant:land_binary_mask', 0)
        o.set_config('general:use_auto_landmask', False)
        o.add_reader(norkyst)
        assert o.mode == Mode.Config
        o.seed_elements(lon=3, lat=60, time=datetime.now())
        assert o.mode == Mode.Ready
        with self.assertRaises(WrongMode):  # Cannot add readers after elements have been seeded
            o.add_reader(norkyst)
        with self.assertRaises(WrongMode):  # Cannot set config after readers have been added
            o.set_config('seed:ocean_only', False)
        o.run(steps=1)
        assert o.mode == Mode.Result


    def test_seed_cone(self):
        o = OceanDrift(loglevel=20)
        o.seed_cone(time=[datetime.now(),
                datetime.now() + timedelta(hours=3)],
                number=100, lat=[60.5, 60.6], lon=[4.4, 4.5])
        self.assertAlmostEqual(o.elements_scheduled.lon[50], 4.450, 2)

    def test_seed_outside_coverage(self):
        """Test seeding"""
        o = OpenOil(loglevel=0)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        landmask = reader_global_landmask.Reader()
        o.add_reader([landmask, norkyst])
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('seed:oil_type', 'SNORRE B 2004')
        o.seed_elements(5, 63, number=5,
                        time=norkyst.start_time - 24*timedelta(hours=24))
        # Check that the oiltype is taken from config
        self.assertEqual(o.oil_name, o.get_config('seed:oil_type'))
        self.assertEqual(o.oil_name, 'SNORRE B 2004')
        with self.assertRaises(ValueError):
            o.run(steps=3, time_step=timedelta(minutes=15))

    def test_invalid_config(self):
        o = OceanDrift(loglevel=20)
        with self.assertRaises(ValueError):  # outside min-max
            o.set_config('seed:number', 0)
        with self.assertRaises(ValueError):  # not in list/enum
            o.set_config('vertical_mixing:diffusivitymodel', 'not_in_list')
        o.set_config('seed:number', 100)
        self.assertEqual(o.get_config('seed:number'), 100)

    def test_config_suggestion(self):
        o = Leeway(loglevel=20)
        try:
            o.set_config('seed:object_type', 'person')
        except Exception as e:
            self.assertTrue('Did you mean' in str(e))

    def test_config_seed(self):
        #o = Leeway(loglevel=20)
        o = OceanDrift(loglevel=20)
        o.list_configspec()

    def test_config_constant_fallback(self):
        o = OceanDrift(loglevel=0)
        reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
        reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        print(reader_norkyst, reader_arome)
        o.add_reader([reader_norkyst, reader_arome])
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 1)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.set_config('environment:constant:x_wind', 0)
        o.set_config('environment:constant:y_wind', 5)
        o.seed_elements(lon=4, lat=60, time=reader_arome.end_time -
                        timedelta(hours=3), number=1)
        o.run(duration=timedelta(hours=6))
        y_wind = np.array(o.get_property('y_wind')[0][:,0])
        x_current = np.array(o.get_property('x_sea_water_velocity')[0][:,0])
        # Check that constant wind is used for whole simulation
        self.assertAlmostEqual(y_wind[0], 5, 2)
        self.assertAlmostEqual(y_wind[-1], 5, 2)
        # Check that fallback current is used only after end of reader
        self.assertAlmostEqual(x_current[0], 0.155, 2)
        self.assertAlmostEqual(x_current[-1], 1, 2)

    def test_runge_kutta(self):
        number = 50
        # With Euler
        o = OceanDrift(loglevel=30, seed=0)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.add_reader([norkyst])
        z=-40*np.random.rand(number)
        o.seed_elements(5, 62.5, number=number, radius=5000, z=z,
                        time=norkyst.start_time)
        o.run(steps=4*3, time_step=timedelta(minutes=15))
        # With Runge-Kutta
        o2 = OceanDrift(loglevel=30, seed=0)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o2.set_config('environment:fallback:land_binary_mask', 0)
        o2.add_reader([norkyst])
        o2.set_config('drift:advection_scheme', 'runge-kutta')
        z=-40*np.random.rand(number)
        o2.seed_elements(5, 62.5, number=number, radius=5000, z=z,
                        time=norkyst.start_time)
        o2.run(steps=4*3, time_step=timedelta(minutes=15))
        # And finally repeating the initial run to check that indetical
        o3 = OceanDrift(loglevel=30, seed=0)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o3.set_config('environment:fallback:land_binary_mask', 0)
        o3.add_reader([norkyst])
        z=-40*np.random.rand(number)
        o3.seed_elements(5, 62.5, number=number, radius=5000, z=z,
                        time=norkyst.start_time)
        o3.run(steps=4*3, time_step=timedelta(minutes=15))
        # Check that we get some difference with Runge-Kutta:
        self.assertIsNone(np.testing.assert_array_almost_equal(
            (o2.elements.lon-o.elements.lon).max(), 0.0015, 3))
            #(o2.elements.lon-o.elements.lon).max(), 0.013, 3))
        # Check that runs with Euler are identical
        self.assertIsNone(np.testing.assert_array_almost_equal(
            (o3.elements.lon-o.elements.lon).max(), 0))

    def test_seed_polygon(self):
        o = OpenOil(loglevel=0)
        number = 10
        lonvec = np.array([2, 3, 3, 2])
        latvec = np.array([60, 60, 61, 61])
        time=datetime(2015, 1, 1, 12, 5, 17)
        o.set_config('seed:oil_type', 'HEIDRUN')
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.seed_within_polygon(lonvec, latvec, number=number,
                              time=time, wind_drift_factor=.09)
        self.assertEqual(o.num_elements_scheduled(), number)
        self.assertEqual(o.elements_scheduled_time[0], time)
        self.assertAlmostEqual(o.elements_scheduled.wind_drift_factor, .09)
        # Check that oil type is taken fom config
        self.assertEqual(o.oil_name, 'HEIDRUN')

    def test_seed_polygon_timespan(self):
        o = OceanDrift(loglevel=20)
        number = 10
        lonvec = np.array([2, 3, 3, 2])
        latvec = np.array([60, 60, 61, 61])
        time=[datetime(2015, 1, 1, 12, 5, 17),
              datetime(2015, 1, 1, 18, 5, 17)]
        o.seed_within_polygon(lonvec, latvec, number=number, time=time)
        self.assertEqual(o.num_elements_scheduled(), number)
        self.assertEqual(o.elements_scheduled_time[0], time[0])
        self.assertEqual(o.elements_scheduled_time[-1], time[-1])

    @unittest.skipIf(has_ogr is False,
                     'OGR library needed to parse WKT')
    def test_seed_wkt(self):
        wkt = 'MULTIPOLYGON(((7.784889 64.353442,7.777561 64.353842,7.774236 64.354707,7.770215 64.355829,7.774269 64.356015,7.776829 64.356863,7.779107 64.3578,7.782827 64.358355,7.786346 64.359615,7.787109 64.361975,7.790125 64.361132,7.794584 64.359908,7.798455 64.359624,7.797258 64.358193,7.79978 64.356904,7.795957 64.356494,7.792955 64.355335,7.789134 64.355339,7.784889 64.353442)))'
        o = OceanDrift(loglevel=20)
        o.set_config('seed:number', 100)
        o.seed_from_wkt(wkt, time=datetime.now())
        wkt_multi = 'MULTIPOLYGON(((2.458058 59.178919,2.456276 59.179283,2.454867 59.180692,2.45277 59.182852,2.452521 59.183759,2.452675 59.184726,2.451365 59.18534,2.451436 59.186609,2.450835 59.188138,2.449576 59.189435,2.447393 59.190818,2.447211 59.191915,2.446273 59.193573,2.445551 59.19423,2.446597 59.195015,2.44838 59.194651,2.450277 59.193,2.452377 59.191919,2.453315 59.19026,2.45457 59.187885,2.455473 59.186131,2.457033 59.18461,2.458774 59.181992,2.458971 59.180403,2.459775 59.179444,2.459606 59.178969,2.458058 59.178919)),((2.442682 59.197444,2.440531 59.198922,2.439575 59.199994,2.440874 59.200951,2.439596 59.20166,2.436232 59.202958,2.433255 59.203728,2.42982 59.203756,2.428 59.202946,2.425857 59.200693,2.42454 59.199149,2.422418 59.198563,2.419404 59.198158,2.417332 59.197175,2.41514 59.19532,2.412395 59.194596,2.410072 59.194519,2.409481 59.193397,2.408199 59.191947,2.405959 59.190489,2.403129 59.188988,2.401292 59.18759,2.398331 59.187867,2.395639 59.187825,2.393585 59.187428,2.389665 59.187697,2.38736 59.188208,2.386923 59.189132,2.390625 59.188785,2.392191 59.189424,2.395825 59.188887,2.398602 59.188627,2.402104 59.189869,2.403773 59.191871,2.407276 59.193113,2.407648 59.194158,2.407751 59.195522,2.410008 59.196488,2.411979 59.197187,2.41439 59.19912,2.415839 59.199965,2.417946 59.201043,2.417796 59.202235,2.414886 59.203195,2.411923 59.203473,2.40923 59.203431,2.409753 59.204363,2.412549 59.20469,2.415342 59.203937,2.41891 59.20321,2.420325 59.203961,2.420463 59.20542,2.419357 59.207683,2.4218 59.208631,2.420303 59.209262,2.418925 59.210766,2.421401 59.21073,2.424984 59.20951,2.425201 59.208508,2.425939 59.207359,2.428832 59.205812,2.431004 59.206001,2.433124 59.205507,2.436926 59.204365,2.439568 59.203724,2.441518 59.202755,2.442879 59.201744,2.443246 59.20063,2.443311 59.199741,2.444589 59.199032,2.445428 59.198168,2.445088 59.197218,2.442682 59.197444)))'
        o.seed_from_wkt(wkt_multi, time=datetime.now(), number=200)
        self.assertEqual(len(o.elements_scheduled), 300)
        self.assertAlmostEqual(
            o.elements_scheduled.lat.max(), 64.36, 2)
        self.assertAlmostEqual(
            o.elements_scheduled.lat.min(), 59.18, 2)

    @unittest.skipIf(has_ogr is False,
                     'OGR library needed to read shapefiles')
    def test_seed_shapefile(self):
        o = OceanDrift(loglevel=20)
        o.seed_from_shapefile(o.test_data_folder() +
                                  'shapefile_spawning_areas/Torsk.shp',
                                  number=100, layername=None,
                                  featurenum=[2, 4], time=datetime.now())
        self.assertEqual(len(o.elements_scheduled), 100)
        o.seed_from_shapefile(o.test_data_folder() +
                                  'shapefile_spawning_areas/Torsk.shp',
                                  number=300, layername=None,
                                  featurenum=None, time=datetime.now())
        self.assertEqual(len(o.elements_scheduled), 400)
        self.assertAlmostEqual(o.elements_scheduled.lat[-1], 52.5, 2)

    #@unittest.skipIf(has_ogr is False,
    #                 'GDAL library needed to read shapefiles')
    #def test_write_geotiff(self):
    #    o = OceanDrift(loglevel=20)
    #    o.seed_elements(lon=4, lat=60, time=datetime(2016, 1, 1))
    #    o.run(steps=3)
    #    o.write_geotiff('geotiff.tif')

    def test_seed_single_point_over_time(self):
        """Test a model run"""
        self.make_OceanDrift_object()
        self.o.seed_elements(2.0, 61.0, radius=0, number=9,
                             time=[datetime(2015, 1, 1), datetime(2015, 1, 3)])

        # Check that 6 elements are scheduled, but none seeded
        self.assertEqual(self.o.num_elements_scheduled(), 9)
        self.assertEqual(self.o.num_elements_active(), 0)
        self.assertEqual(self.o.num_elements_activated(), 0)
        self.assertEqual(self.o.num_elements_deactivated(), 0)
        self.assertEqual(self.o.num_elements_total(), 9)
        # Run simulation
        self.o.run(steps=30)
        # Check that 1 element is deactivated (stranded),
        # 1 yet not seeded and 4 active
        self.assertEqual(self.o.num_elements_scheduled(), 4)
        self.assertEqual(self.o.num_elements_active(), 5)
        self.assertEqual(self.o.num_elements_activated(), 5)
        self.assertEqual(self.o.num_elements_deactivated(), 0)
        self.assertEqual(self.o.num_elements_total(), 9)

    def test_temporal_seed(self):
        self.o = OceanDrift(loglevel=20)
        self.o.set_config('environment:fallback:x_sea_water_velocity', 1)
        self.o.set_config('environment:fallback:land_binary_mask', 0)
        # Seed elements on a grid at regular time interval
        start_time = datetime(2016,9,16)
        time_step = timedelta(hours=6)
        num_steps = 10
        lon = 4.4
        lat = 60.0
        for i in range(num_steps+1):
            self.o.seed_elements(lon, lat, radius=0, number=2,
                                 time=start_time + i*time_step)
        # Running model
        self.o.run(steps=20, time_step=3600, outfile='temporal_seed.nc')
        self.o = OceanDrift(loglevel=20)
        # Check that data imported is properly masked
        self.o.io_import_file('temporal_seed.nc')
        self.assertTrue(self.o.history['lon'].max() < 1000)
        self.assertTrue(self.o.history['lon'].min() > -1000)
        self.assertTrue(self.o.history['lon'].mask[5,5])
        self.assertFalse(self.o.history['lon'].mask[1,1])
        os.remove('temporal_seed.nc')

    def test_vertical_mixing(self):
        # Export to file only at end
        o1 = PelagicEggDrift(loglevel=20)  # Profiles and vertical mixing
        norkyst = reader_netCDF_CF_generic.Reader(o1.test_data_folder() +
          '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o1.add_reader([norkyst])
        o1.set_config('environment:fallback:x_wind', 8)
        o1.set_config('environment:fallback:land_binary_mask', 0)
        o1.set_config('vertical_mixing:timestep', 20.) # seconds
        o1.seed_elements(4.1, 63.3, radius=1000, number=100,
                         time=norkyst.start_time)

        o1.run(steps=20, time_step=300, time_step_output=1800,
               export_buffer_length=10, outfile='verticalmixing.nc')
        self.assertAlmostEqual(o1.history['z'].min(), -43.67, 1)
        self.assertAlmostEqual(o1.history['z'].max(), 0.0, 1)
        os.remove('verticalmixing.nc')

    def test_vertical_mixing_profiles(self):
        # Testing an isolated mixing timestep

        cases = [  # Some cases with expected outcome
                {'vt': 0, 'K': 0, 'K_below': .01, 'T': 60, # No mixing
                    'zmin': -10, 'zmax': -10, 'zmean': -10},
                {'vt': -.005, 'K': 0, 'K_below': .01, 'T': 60, # Sinking
                    #'zmin': -74.0, 'zmax': -21.4, 'zmean': -50.2},  # With old seed_elements
                    'zmin': -74.79, 'zmax': -21.6, 'zmean': -49.97},
                {'vt': 0, 'K': .01, 'K_below': .01, 'T': 60, # Mixing
                    #'zmin': -39.8, 'zmax': -0.1, 'zmean': -14.5},
                    'zmin': -42.76, 'zmax': -0.02, 'zmean': -14.38},
                {'vt': .005, 'K': .01, 'K_below': .01, 'T': 60, # Mixing and rising
                    #'zmin': -8.1, 'zmax': -0.01, 'zmean': -2.1},
                    'zmin': -7.86, 'zmax': -0.01, 'zmean': -2.1},
                {'vt': -0.005, 'K': .01, 'K_below': .01, 'T': 60, # Mixing and sinking
                    #'zmin': -75.8, 'zmax': -20.7, 'zmean': -48.1},
                    'zmin': -78.76, 'zmax': -19.74, 'zmean': -48.0},
                {'vt': 0, 'K': .02, 'K_below': .001, 'T': 60,  # Mixing in mixed layer
                    #'zmin': -22.8, 'zmax': -0.1, 'zmean': -9.8},
                    'zmin': -21.3, 'zmax': -0.1, 'zmean': -9.55},
                ]

        N=100
        z = np.arange(0, -30, -2)
        time = datetime.now()
        for case in cases:
            diffusivity = np.ones(z.shape)*case['K']
            diffusivity[z<-15] = case['K_below']
            o = OceanDrift(loglevel=20)
            o.set_config('drift:vertical_mixing', True)
            o.set_config('vertical_mixing:diffusivitymodel', 'environment')
            o.set_config('vertical_mixing:timestep', case['T'])
            o.set_config('environment:fallback:land_binary_mask', 0)
            o.seed_elements(lon=4, lat=60, z=-10, time=time, number=N,
                            terminal_velocity=case['vt'])
            o.time = time
            o.time_step = timedelta(hours=2)
            o.release_elements()
            o.environment = np.array([(100, 0) for _ in range(N)],
                         dtype=[('sea_floor_depth_below_sea_level', np.float32),
                                ('sea_surface_height', np.float32)]).view(np.recarray)
            o.environment.ocean_mixed_layer_thickness = np.ones(N)*50
            o.environment_profiles = {'z': z, 'ocean_vertical_diffusivity':
                                      np.tile(diffusivity, (N, 1)).T}
            o.env.finalize()
            o.vertical_mixing()
            self.assertAlmostEqual(o.elements.z.min(), case['zmin'], 1)
            self.assertAlmostEqual(o.elements.z.max(), case['zmax'], 1)
            self.assertAlmostEqual(o.elements.z.mean(), case['zmean'], 1)

    def test_export_step_interval(self):
        # Export to file only at end
        o1 = OceanDrift(loglevel=20)
        norkyst = reader_netCDF_CF_generic.Reader(o1.test_data_folder() +
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        o1.add_reader(norkyst)
        o1.set_config('environment:fallback:land_binary_mask', 0)
        o1.seed_elements(4.25, 60.2, radius=1000, number=10,
                        time=norkyst.start_time)
        o1.run(steps=40)
        # Export to file during simulation
        o2 = OceanDrift(loglevel=20)
        o2.add_reader(norkyst)
        o2.set_config('environment:fallback:land_binary_mask', 0)
        o2.seed_elements(4.25, 60.2, radius=1000, number=10,
                        time=norkyst.start_time)
        o2.run(steps=40, export_buffer_length=6,
               outfile='export_step_interval.nc')
        self.assertIsNone(np.testing.assert_array_equal(
            o1.history['lon'].compressed(),
            o2.history['lon'].compressed()))
        # Finally check when steps is multiple of export_buffer_length
        o3 = OceanDrift(loglevel=20)
        o3.add_reader(norkyst)
        o3.set_config('environment:fallback:land_binary_mask', 0)
        o3.seed_elements(4.25, 60.2, radius=1000, number=10,
                        time=norkyst.start_time)
        o3.run(steps=42)
        # Export to file during simulation
        o4 = OceanDrift(loglevel=20)
        o4.add_reader(norkyst)
        o4.set_config('environment:fallback:land_binary_mask', 0)
        o4.seed_elements(4.25, 60.2, radius=1000, number=10,
                        time=norkyst.start_time)
        o4.run(steps=42, export_buffer_length=6,
               outfile='export_step_interval.nc')
        self.assertIsNone(np.testing.assert_array_equal(
            o3.history['lon'].compressed(),
            o4.history['lon'].compressed()))
        os.remove('export_step_interval.nc')

    def test_export_final_timestep(self):
        o = OceanDrift()
        o.set_config('environment:constant:land_binary_mask', 0)
        o.set_config('general:use_auto_landmask', False)
        o.seed_elements(lon=0, lat=0, radius=500, number=100,
                        time=[datetime(2010,1,1), datetime(2010,1,3)])
        o.run(duration=timedelta(hours=20), time_step=3600, time_step_output=3600*3)
        index_of_first, index_of_last = \
            o.index_of_activation_and_deactivation()
        assert o.num_elements_active() == len(index_of_first)

    def test_buffer_length_stranding(self):
        o1 = OceanDrift(loglevel=30)
        norkyst = reader_netCDF_CF_generic.Reader(o1.test_data_folder() +
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        landmask = reader_global_landmask.Reader()
        o1.add_reader([landmask])
        o1.set_config('environment:fallback:x_sea_water_velocity', 0.8)  # onshore drift
        o1.seed_elements(4.8, 60.2, radius=5000, number=100,
                        time=norkyst.start_time)
        o1.run(steps=100,
               time_step=900,
               time_step_output=3600,
               export_buffer_length=10)
        # Without buffer
        o2 = OceanDrift(loglevel=30)
        o2.add_reader([landmask])
        o2.set_config('environment:fallback:x_sea_water_velocity', 0.8)  # onshore drift
        o2.seed_elements(4.8, 60.2, radius=5000, number=100,
                        time=norkyst.start_time)
        o2.run(steps=100,
               time_step=900,
               time_step_output=3600,
               outfile='test_buffer_length_stranding.nc')
        self.assertIsNone(np.testing.assert_array_equal(
            o1.history['lon'].compressed(),
            o2.history['lon'].compressed()))
        self.assertIsNone(np.testing.assert_array_almost_equal(
            o1.history['status'].compressed(),
            o2.history['status'].compressed()))
        os.remove('test_buffer_length_stranding.nc')

    def test_output_time_step(self):
        o1 = OceanDrift(loglevel=30)
        norkyst = reader_netCDF_CF_generic.Reader(o1.test_data_folder() +
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        landmask = reader_global_landmask.Reader()
        o1.add_reader([landmask, norkyst])
        o1.seed_elements(4.96, 60.1, radius=3000, number=100,
                        time=norkyst.start_time)
        o1.run(duration=timedelta(hours=12),
                   time_step=timedelta(minutes=30),
                   time_step_output=timedelta(minutes=30),
                   outfile='test_time_step30.nc')
        # Check length of time array and output array
        time = o1.get_time_array()[0]
        self.assertEqual(o1.history.shape[1], len(time))
        self.assertEqual(o1.start_time, time[0])
        self.assertEqual(o1.time, time[-1])
        # Second run, with larger output time step
        o2 = OceanDrift(loglevel=30)
        o2.add_reader([landmask, norkyst])
        o2.seed_elements(4.96, 60.1, radius=3000, number=100,
                        time=norkyst.start_time)
        o2.run(duration=timedelta(hours=12),
                   time_step=timedelta(minutes=30),
                   time_step_output=timedelta(minutes=60),
                   outfile='test_time_step60.nc')
        self.assertEqual(o1.history.shape, (100,25))
        self.assertEqual(o2.history.shape, (100,13))
        # Check that start and end conditions (longitudes) are idential
        self.assertIsNone(np.testing.assert_array_equal(
            o1.history['lon'][:,24].compressed(),
            o2.history['lon'][:,12].compressed()))
        self.assertIsNone(np.testing.assert_array_equal(
            o1.history['lon'][:,0].compressed(),
            o2.history['lon'][:,0].compressed()))
        # Check that also run imported from file is identical
        o1i = OceanDrift(loglevel=20)
        o1i.io_import_file('test_time_step30.nc')
        o2i = OceanDrift(loglevel=20)
        o2i.io_import_file('test_time_step60.nc')
        os.remove('test_time_step30.nc')
        os.remove('test_time_step60.nc')
        self.assertIsNone(np.testing.assert_array_equal(
            o2i.history['lon'][:,12].compressed(),
            o2.history['lon'][:,12].compressed()))
        # Check number of activated elements
        self.assertEqual(o1.num_elements_total(), o2.num_elements_total())
        self.assertEqual(o1.num_elements_total(), o1i.num_elements_total())
        self.assertEqual(o1.num_elements_total(), o2i.num_elements_total())
        # Check number of deactivated elements
        self.assertEqual(o1.num_elements_deactivated(),
                         o2.num_elements_deactivated())
        self.assertEqual(o1.num_elements_deactivated(),
                         o1i.num_elements_deactivated())
        self.assertEqual(o1.num_elements_deactivated(),
                         o2i.num_elements_deactivated())

    def test_time_step_config(self):
        # Default
        o = OceanDrift(loglevel=50)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.seed_elements(lon=4, lat=60, time=datetime.now())
        o.run(steps=2)
        self.assertEqual(o.time_step.total_seconds(), 3600)
        self.assertEqual(o.time_step_output.total_seconds(), 3600)
        # Setting time_step
        o = OceanDrift(loglevel=50)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.seed_elements(lon=4, lat=60, time=datetime.now())
        o.run(steps=2, time_step=1800)
        self.assertEqual(o.time_step.total_seconds(), 1800)
        # Setting time_step and time_step_output
        o = OceanDrift(loglevel=50)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.seed_elements(lon=4, lat=60, time=datetime.now())
        o.run(steps=2, time_step=1800, time_step_output=3600)
        self.assertEqual(o.time_step.total_seconds(), 1800)
        self.assertEqual(o.time_step_output.total_seconds(), 3600)
        # time_step from config
        o = OceanDrift(loglevel=50)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('general:time_step_minutes', 15)
        o.seed_elements(lon=4, lat=60, time=datetime.now())
        o.run(steps=2)
        self.assertEqual(o.time_step.total_seconds(), 900)
        self.assertEqual(o.time_step_output.total_seconds(), 900)
        # time_step and time_step_output from config
        o = OceanDrift(loglevel=50)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('general:time_step_minutes', 15)
        o.set_config('general:time_step_output_minutes', 120)
        o.seed_elements(lon=4, lat=60, time=datetime.now())
        o.run(steps=2)
        self.assertEqual(o.time_step.total_seconds(), 900)
        self.assertEqual(o.time_step_output.total_seconds(), 7200)


    def test_reader_boundary(self):
        # Check that the element outside reader coverage is
        # not deactivated if fallback value exist
        o = OceanDrift()
        nordic3d = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        lon = [12.0, 12.0]
        lat = [70.0, 70.5]
        o.add_reader(nordic3d)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.seed_elements(lon, lat, number=2, radius=0,
                        time=nordic3d.start_time)
        o.run(steps=2, time_step=3600)
        self.assertEqual(o.num_elements_active(), 2)
        self.assertEqual(o.num_elements_deactivated(), 0)
        # Check that the outside element is deactivated,
        # if no fallback value exists
        o = OceanDrift()
        o.set_config('environment:fallback:x_sea_water_velocity', None)
        o.add_reader(nordic3d)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.seed_elements(lon, lat, number=2, radius=0,
                        time=nordic3d.start_time)
        o.run(steps=2, time_step=3600)
        self.assertEqual(o.num_elements_active(), 1)
        self.assertEqual(o.num_elements_deactivated(), 1)

    def test_seed_seafloor(self):
        o = OpenOil(loglevel=50)
        reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        # Adding reader as lazy, to test seafloor seeding
        o.add_readers_from_list([o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc'])
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)
        #o.set_config('environment:fallback:x_sea_water_velocity', 0)
        #o.set_config('environment:fallback:y_sea_water_velocity', 0)
        lon = 4.5; lat = 62.0
        o.set_config('seed:droplet_diameter_min_subsea', 0.0010)  # s
        o.set_config('seed:droplet_diameter_max_subsea', 0.0010)  # s
        o.set_config('drift:vertical_mixing', True)
        o.set_config('vertical_mixing:timestep', 1)  # s
        o.seed_elements(lon, lat, z='seafloor', time=reader_norkyst.start_time,
                        density=1000, oiltype='GENERIC BUNKER C')
        o.run(steps=3, time_step=300, time_step_output=300)
        #o.plot_property('z')
        z, status = o.get_property('z')
        self.assertAlmostEqual(z[0,0], -147.3, 1)  # Seeded at seafloor depth
        self.assertAlmostEqual(z[-1,0], -132.31, 1)  # After some rising

    def test_seed_above_seafloor(self):
        o = OpenOil(loglevel=30)
        reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.add_reader([reader_norkyst])
        lon = 4.5; lat = 62.0
        o.set_config('seed:droplet_diameter_min_subsea', 0.0010)  # s
        o.set_config('seed:droplet_diameter_max_subsea', 0.0010)  # s
        o.set_config('drift:vertical_mixing', True)
        o.set_config('vertical_mixing:timestep', 1)  # s
        # Seed elements 50 meters above seafloor:
        o.seed_elements(lon, lat, z='seafloor+50', time=reader_norkyst.start_time,
                        density=1000, oil_type='AASGARD A 2003')
        o.run(steps=3, time_step=300, time_step_output=300)
        #o.plot_property('z')
        z, status = o.get_property('z')
        self.assertAlmostEqual(z[0,0], -97.3, 1)  # Seeded at seafloor depth
        self.assertAlmostEqual(z[-1,0], -38.3, 1)  # After some rising

    def test_seed_below_reader_coverage(self):
        o = OpenOil(loglevel=20)
        reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)
        o.add_reader([reader_norkyst])
        lon = 5.0; lat = 64.0
        o.set_config('seed:droplet_diameter_min_subsea', 0.0005)
        o.set_config('seed:droplet_diameter_max_subsea', 0.005)
        #o.set_config('vertical_mixing:TSprofiles', True)
        o.set_config('drift:vertical_mixing', True)
        o.set_config('vertical_mixing:timestep', 1)  # s
        o.seed_elements(lon, lat, z=-350, time=reader_norkyst.start_time,
                        density=1000, oil_type='AASGARD A 2003')
        o.run(steps=3, time_step=300, time_step_output=300)
        z, status = o.get_property('z')
        self.assertAlmostEqual(z[-1,0], -235.5, 1)  # After some rising

    def test_seed_below_seafloor(self):
        o = OpenOil(loglevel=20)
        reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o.add_reader([reader_norkyst])
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        lon = 4.5; lat = 62.0
        o.set_config('seed:droplet_diameter_min_subsea', 0.0005)
        o.set_config('seed:droplet_diameter_max_subsea', 0.001)
        o.set_config('drift:vertical_mixing', True)
        o.set_config('vertical_mixing:timestep', 1)  # s
        o.seed_elements(lon, lat, z=-5000, time=reader_norkyst.start_time,
                        density=1000, oil_type='GENERIC BUNKER C')
        o.run(steps=3, time_step=300, time_step_output=300)
        z, status = o.get_property('z')
        self.assertAlmostEqual(z[0,0], -147.3, 1)  # Seeded at seafloor depth
        self.assertAlmostEqual(z[-1,0], -138.3, 1)  # After some rising

    def test_seed_below_seafloor_deactivating(self):
        o = OpenOil(loglevel=50)
        reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o.add_reader([reader_norkyst])
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        lon = 4.5; lat = 62.0
        o.set_config('seed:droplet_diameter_min_subsea', 0.0005)
        o.set_config('seed:droplet_diameter_max_subsea', 0.001)
        o.set_config('general:seafloor_action', 'deactivate')  # This time we deactivate
        o.set_config('drift:vertical_mixing', True)
        o.set_config('vertical_mixing:timestep', 1)  # s
        o.seed_elements(lon, lat, z=[-5000, -100], time=reader_norkyst.start_time,
                        density=1000, number=2, oil_type='AASGARD A 2003')
        o.run(steps=3, time_step=300, time_step_output=300)
        z, status = o.get_property('z')
        self.assertEqual(o.num_elements_total(), 2)
        self.assertEqual(o.num_elements_active(), 1)
        self.assertEqual(o.num_elements_deactivated(), 1)
        self.assertAlmostEqual(z[0,1], -100, 1)  # Seeded at seafloor depth
        self.assertAlmostEqual(z[-1,1], -56.7, 1)  # After some rising

    def test_lift_above_seafloor(self):
        # See an element at some depth, and progapate towards coast
        # (shallower water) and check that it is not penetrating seafloor
        o = OceanDrift(loglevel=50)
        o.set_config('drift:max_speed', 100)
        reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        reader_norkyst.buffer = 200
        o.add_reader([reader_norkyst],
                     variables='sea_floor_depth_below_sea_level')
        o.set_config('environment:fallback:x_sea_water_velocity', 10) # Pure eastward motion
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('drift:vertical_mixing', False)
        o.seed_elements(3.9, 62.0, z=-200, time=reader_norkyst.start_time)
        o.run(steps=12, time_step=300)
        seafloor_depth, status = o.get_property('sea_floor_depth_below_sea_level')
        z, status = o.get_property('z')

        # Uncomment to plot
        #import matplotlib.pyplot as plt
        #plt.plot(-seafloor_depth, label='Seafloor depth')
        #plt.plot(z, label='Element depth')
        #plt.legend(loc='best')
        #plt.show()

        # Check that element has not penetrated seafloor
        self.assertFalse(o.elements.z <
                         -o.environment.sea_floor_depth_below_sea_level)
        self.assertIsNone(np.testing.assert_array_almost_equal(
            o.elements.z, -134.5, 1))

    def test_seed_on_land(self):
        o = OceanDrift(loglevel=50)
        o.seed_elements(lon=9, lat=60, time=datetime.now(), number=100)
        outfile='out.nc'
        with self.assertRaises(ValueError):
            o.run(steps=4, time_step=1800, time_step_output=3600,
                  outfile=outfile)
        os.remove(outfile)
        #o.write_netcdf_density_map(outfile)
        #os.remove(outfile)

    def test_retirement(self):
        o = OceanDrift(loglevel=0)
        o.set_config('drift:max_age_seconds', 5000)
        o.set_config('environment:fallback:x_sea_water_velocity', .5)
        o.set_config('environment:fallback:y_sea_water_velocity', .3)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.seed_elements(lon=0, lat=60, number=10,
                        time=[datetime.now(),
                              datetime.now() + timedelta(seconds=6000)])

        o.run(time_step=1000, duration=timedelta(seconds=7000))
        self.assertEqual(o.num_elements_deactivated(), 5)

    def test_outside_domain(self):
        o = OceanDrift(loglevel=50)
        reader_osc_x = reader_oscillating.Reader(
                'x_sea_water_velocity', amplitude=1,
                zero_time=datetime.now())
        reader_osc_y = reader_oscillating.Reader(
                'y_sea_water_velocity', amplitude=1,
                zero_time=datetime.now())
        o.add_reader([reader_osc_x, reader_osc_y])
        o.set_config('drift:deactivate_east_of', 2.1)
        o.set_config('drift:deactivate_west_of', 1.9)
        o.set_config('drift:deactivate_south_of', 59.9)
        o.set_config('drift:deactivate_north_of', 60.1)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.seed_elements(lon=2, lat=60, number=1000,
                        time=datetime.now(), radius=10000)
        o.run(duration=timedelta(hours=5))
        self.assertEqual(o.num_elements_deactivated(), 768)
        self.assertEqual(o.num_elements_active(), 232)

    def test_seed_time_backwards_run(self):
        o = OceanDrift(loglevel=20)
        o.set_config('drift:max_age_seconds', 2000)
        o.set_config('environment:fallback:x_sea_water_velocity', .5)
        o.set_config('environment:fallback:y_sea_water_velocity', .3)
        o.set_config('environment:fallback:land_binary_mask', 0)
        time = [datetime(2018,1,1,i) for i in range(10)]
        o.seed_elements(lon=0, lat=60, time=time)
        o.seed_elements(lon=1, lat=60, time=datetime(2018,1,1,7))
        o.run(end_time=datetime(2018,1,1,2), time_step=-1800)
        self.assertEqual(o.num_elements_scheduled(), 3)
        self.assertEqual(o.num_elements_active(), 8)
        self.assertEqual(o.steps_calculation, 14)

    def test_oil_mixed_to_seafloor(self):
        o = OpenOil(loglevel=30)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o.add_reader(norkyst)
        o.set_config('processes:evaporation', False)
        o.set_config('environment:fallback:x_wind', 25)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:ocean_vertical_diffusivity', 0.9)
        o.seed_elements(lon=5.38, lat=62.77, time=norkyst.start_time,
                        number=100, radius=5000)
        o.run(end_time=norkyst.end_time)
        self.assertEqual(o.num_elements_active(), 100)

    def test_unseeded_elements(self):
        o = PlastDrift()
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 1)
        # Seeding elements for 12 hours, but running only 6
        time = datetime(2019, 8, 30, 12)
        o.seed_elements(lon=4.85, lat=60, number=10,
                        time=[time, time + timedelta(hours=6)],
                        origin_marker=7)
        o.seed_elements(lon=4.75, lat=60, number=10,
                        time=[time, time + timedelta(hours=6)],
                        origin_marker=8)
        o.run(duration=timedelta(hours=3))
        self.assertEqual(o.history.shape[0], 10)
        self.assertEqual(o.history.shape[1], 4)
        self.assertEqual(o.history['origin_marker'].min(), 7)
        self.assertEqual(o.history['origin_marker'].max(), 8)

    def test_no_active_but_still_unseeded_elements(self):
        o = OceanDrift(loglevel=20)
        # deactivate elements after 3 hours
        o.set_config('drift:max_age_seconds', 3600*3)
        o.set_config('environment:fallback:land_binary_mask', 0)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o.add_reader(norkyst)
        # seed two elements at 6 hour interval
        o.seed_elements(number=2, lon=4, lat=62,
            time=[norkyst.start_time, norkyst.start_time+timedelta(hours=6)])
        o.run(duration=timedelta(hours=8), outfile='test.nc')
        os.remove('test.nc')
        # Check that simulations has run until scheduled end
        self.assertEqual(o.steps_calculation, 8)

@pytest.mark.slow
def test_plot_animation(tmpdir):
    o = OceanDrift(loglevel=0)
    o.set_config('environment:fallback:x_sea_water_velocity', .5)
    o.set_config('environment:fallback:y_sea_water_velocity', .3)
    o.seed_elements(lon=3, lat=60, radius=1000,
                    time=datetime.now(), number=100)
    o.run(steps=5)
    # Plot
    o.plot(filename='%s/test_plot.png' % tmpdir, lscale='c')
    assert os.path.exists('%s/test_plot.png' % tmpdir)
    # Animation
    # Temporarily skipping mp4 due to bug in Conda
    #o.animation(filename='test_plot.mp4')
    #assert os.path.exists('test_plot.mp4')
    #os.remove('test_plot.mp4')


if __name__ == '__main__':
    unittest.main()
