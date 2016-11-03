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
import os
import inspect

import numpy as np

from opendrift.readers import reader_ArtificialOceanEddy
from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.oceandrift3D import OceanDrift3D
from opendrift.models.openoil3D import OpenOil3D
from opendrift.models.pelagicegg import PelagicEggDrift

try:
    import ogr
    import osr
    import gdal
    has_ogr = True
except Exception as e:
    print 'GDAL is not available:'
    print e
    has_ogr = False

class TestRun(unittest.TestCase):
    """Tests for (non-scalar) LagrangianArray"""

    def make_OceanDrift_object(self):
        self.o = OceanDrift(loglevel=20)
        self.fake_eddy = reader_ArtificialOceanEddy.Reader(2, 62)
        self.o.use_block = False
        self.o.config['runge_kutta'] = False
        self.reader_basemap = reader_basemap_landmask.Reader(
            llcrnrlon=-1.5, llcrnrlat=59,
            urcrnrlon=7, urcrnrlat=64, resolution='i')
        self.o.add_reader([self.fake_eddy, self.reader_basemap])

    def test_config(self):
        """Test invalid configuration"""
        o = OceanDrift(loglevel=20)
        o.config['drift']['max_age_seconds'] = -1  # invalid value
        o.fallback_values['land_binary_mask'] = 0
        o.seed_elements(5, 60, time=datetime(2016,1,1))
        self.assertRaises(ValueError, o.run, steps=3)

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

    def test_seed_outside_coverage(self):
        """Test seeding"""
        o = OpenOil3D(loglevel=0)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        basemap = reader_basemap_landmask.Reader(
            llcrnrlon=4, llcrnrlat=60, urcrnrlon=6, urcrnrlat=64,
            resolution='c', projection='merc')
        o.add_reader([basemap, norkyst])
        o.seed_elements(5, 63, number=5,
                        time=norkyst.start_time - 24*timedelta(hours=24))
        o.run(steps=3, time_step=timedelta(minutes=15))

    def test_runge_kutta(self):
        number = 50
        # With Euler
        o = OceanDrift3D(loglevel=0, seed=0)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o.fallback_values['land_binary_mask'] = 0
        o.add_reader([norkyst])
        z=-40*np.random.rand(number)
        o.seed_elements(5, 62.5, number=number, radius=5000, z=z,
                        time=norkyst.start_time)
        o.run(steps=4*3, time_step=timedelta(minutes=15))
        # With Runge-Kutta
        o2 = OceanDrift3D(loglevel=20, seed=0)
        o2.fallback_values['land_binary_mask'] = 0
        o2.add_reader([norkyst])
        z=-40*np.random.rand(number)
        o2.seed_elements(5, 62.5, number=number, radius=5000, z=z,
                        time=norkyst.start_time)
        o2.config['drift']['scheme'] = 'runge-kutta'
        o2.run(steps=4*3, time_step=timedelta(minutes=15))
        # And finally repeating the initial run to check that indetical
        o3 = OceanDrift3D(loglevel=20, seed=0)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o3.fallback_values['land_binary_mask'] = 0
        o3.add_reader([norkyst])
        z=-40*np.random.rand(number)
        o3.seed_elements(5, 62.5, number=number, radius=5000, z=z,
                        time=norkyst.start_time)
        o3.run(steps=4*3, time_step=timedelta(minutes=15))
        # Check that we get some difference with Runge-Kutta:
        self.assertAlmostEqual((o2.elements.lon-o.elements.lon).max(),
                                0.0015, 3)
        # Check that runs with Euler are identical
        self.assertEqual((o3.elements.lon-o.elements.lon).max(), 0)

    def test_seed_polygon(self):
        o = OceanDrift(loglevel=20)
        number = 10
        lonvec = np.array([2, 3, 3, 2])
        latvec = np.array([60, 60, 61, 61])
        time=datetime(2015, 1, 1, 12, 5, 17)
        o.seed_within_polygon(lonvec, latvec, number=number, time=time,
                              wind_drift_factor=.09)
        self.assertEqual(o.num_elements_scheduled(), number)
        self.assertEqual(o.elements_scheduled_time[0], time)
        self.assertAlmostEqual(o.elements_scheduled.wind_drift_factor, .09)

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

    #@unittest.skipIf(has_ogr is False,
    #                 'GDAL library needed to read shapefiles')
    #def test_write_geotiff(self):
    #    o = OceanDrift(loglevel=20)
    #    o.seed_elements(lon=4, lat=60, time=datetime(2016, 1, 1))
    #    o.run(steps=3)
    #    o.write_geotiff('geotiff.tif')

    def test1_seed_single_point_over_time(self):
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
        self.o.run(steps=30, outfile='unittest.nc')
        # Check that 1 element is deactivated (stranded),
        # 1 yet not seeded and 4 active
        self.assertEqual(self.o.num_elements_scheduled(), 4)
        self.assertEqual(self.o.num_elements_active(), 5)
        self.assertEqual(self.o.num_elements_activated(), 5)
        self.assertEqual(self.o.num_elements_deactivated(), 0)
        self.assertEqual(self.o.num_elements_total(), 9)

    def test2_seed_elementss(self):
        """Test a model run"""
        self.make_OceanDrift_object()
        self.o.seed_elements([2.0, 4.5, 3.0], [61.0, 60.0, 62.0],
                             radius=0, number=9,
                             time=[datetime(2015, 1, 1), datetime(2015, 1, 3)])

        # Check that 6 elements are scheduled, but none seeded
        self.assertEqual(self.o.num_elements_scheduled(), 9)
        self.assertEqual(self.o.num_elements_active(), 0)
        self.assertEqual(self.o.num_elements_activated(), 0)
        self.assertEqual(self.o.num_elements_deactivated(), 0)
        self.assertEqual(self.o.num_elements_total(), 9)
        # Run simulation
        self.o.run(steps=30, outfile='unittest.nc')
        # Check that 1 element is deactivated (stranded),
        # 1 yet not seeded and 4 active
        self.assertEqual(self.o.num_elements_scheduled(), 4)
        self.assertEqual(self.o.num_elements_active(), 4)
        self.assertEqual(self.o.num_elements_activated(), 5)
        self.assertEqual(self.o.num_elements_deactivated(), 1)
        self.assertEqual(self.o.num_elements_total(), 9)

    def test3_run_import(self):
        """Import output file from previous test, and check elements"""
        self.o = OceanDrift(loglevel=20)
        self.o.io_import_file('unittest.nc')
        self.assertEqual(self.o.num_elements_active(), 4)
        self.assertEqual(self.o.num_elements_activated(), 5)
        self.assertEqual(self.o.num_elements_deactivated(), 1)
        self.assertEqual(self.o.num_elements_total(), 5)

    def test4_cleaning(self):
        """Cleaning up"""
        os.remove('unittest.nc')

    def test_temporal_seed(self):
        self.o = OceanDrift(loglevel=20)
        self.o.fallback_values['x_sea_water_velocity'] = 1
        self.o.fallback_values['land_binary_mask'] = 0
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
        o1 = PelagicEggDrift(loglevel=0)  # Profiles and vertical mixing
        norkyst = reader_netCDF_CF_generic.Reader(o1.test_data_folder() +
          '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        o1.add_reader([norkyst])
        o1.fallback_values['x_wind'] = 8
        o1.fallback_values['land_binary_mask'] = 0
        o1.seed_elements(4.1, 63.3, radius=1000, number=100,
                         time=norkyst.start_time)
        o1.config['turbulentmixing']['timestep'] = 20. # seconds
        o1.run(steps=20, time_step=300, time_step_output=1800,
               export_buffer_length=10, outfile='verticalmixing.nc')
        self.assertAlmostEqual(o1.history['z'].min(), -25.0)
        self.assertAlmostEqual(o1.history['z'].max(), 0.0)
        os.remove('verticalmixing.nc')

    def test_export_step_interval(self):
        # Export to file only at end
        o1 = OceanDrift(loglevel=20)
        norkyst = reader_netCDF_CF_generic.Reader(o1.test_data_folder() +
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        o1.add_reader(norkyst)
        o1.fallback_values['land_binary_mask'] = 0
        o1.seed_elements(4.25, 60.2, radius=1000, number=10,
                        time=norkyst.start_time)
        o1.run(steps=40)
        # Export to file during simulation
        o2 = OceanDrift(loglevel=20)
        o2.add_reader(norkyst)
        o2.fallback_values['land_binary_mask'] = 0
        o2.seed_elements(4.25, 60.2, radius=1000, number=10,
                        time=norkyst.start_time)
        o2.run(steps=40, export_buffer_length=6,
               outfile='export_step_interval.nc')
        self.assertItemsEqual(o1.history['lon'].compressed(),
                              o2.history['lon'].compressed())
        # Finally check when steps is multiple of export_buffer_length
        o3 = OceanDrift(loglevel=20)
        o3.add_reader(norkyst)
        o3.fallback_values['land_binary_mask'] = 0
        o3.seed_elements(4.25, 60.2, radius=1000, number=10,
                        time=norkyst.start_time)
        o3.run(steps=42)
        # Export to file during simulation
        o4 = OceanDrift(loglevel=20)
        o4.add_reader(norkyst)
        o4.fallback_values['land_binary_mask'] = 0
        o4.seed_elements(4.25, 60.2, radius=1000, number=10,
                        time=norkyst.start_time)
        o4.run(steps=42, export_buffer_length=6,
               outfile='export_step_interval.nc')
        self.assertItemsEqual(o3.history['lon'].compressed(),
                              o4.history['lon'].compressed())
        os.remove('export_step_interval.nc')

    def test_buffer_length_stranding(self):
        o1 = OceanDrift(loglevel=0)
        norkyst = reader_netCDF_CF_generic.Reader(o1.test_data_folder() +
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        basemap = reader_basemap_landmask.Reader(
            llcrnrlon=1, llcrnrlat=59.8, urcrnrlon=6, urcrnrlat=61,
            resolution='i', projection='merc')
        o1.add_reader([basemap])
        o1.fallback_values['x_sea_water_velocity'] = 1.0  # onshore drift
        o1.seed_elements(4.8, 60.2, radius=5000, number=100,
                        time=norkyst.start_time)
        o1.run(steps=100,
               time_step=900,
               time_step_output=3600,
               export_buffer_length=10)
        # Without buffer
        o2 = OceanDrift(loglevel=0)
        o2.add_reader([basemap])
        o2.fallback_values['x_sea_water_velocity'] = 1.0  # onshore drift
        o2.seed_elements(4.8, 60.2, radius=5000, number=100,
                        time=norkyst.start_time)
        o2.run(steps=100,
               time_step=900,
               time_step_output=3600,
               outfile='test_buffer_length_stranding.nc')
        self.assertItemsEqual(o1.history['lon'].compressed(),
                              o2.history['lon'].compressed())
        self.assertItemsEqual(o1.history['status'].compressed(),
                              o2.history['status'].compressed())
        os.remove('test_buffer_length_stranding.nc')

    def test_output_time_step(self):
        o1 = OceanDrift(loglevel=0)
        norkyst = reader_netCDF_CF_generic.Reader(o1.test_data_folder() +
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        basemap = reader_basemap_landmask.Reader(
            llcrnrlon=4, llcrnrlat=59.8, urcrnrlon=6, urcrnrlat=61,
            resolution='i', projection='merc')
        o1.add_reader([basemap, norkyst])
        o1.seed_elements(4.95, 60.1, radius=3000, number=100,
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
        o2 = OceanDrift(loglevel=20)
        o2.add_reader([basemap, norkyst])
        o2.seed_elements(4.95, 60.1, radius=3000, number=100,
                        time=norkyst.start_time)
        o2.run(duration=timedelta(hours=12),
                   time_step=timedelta(minutes=30),
                   time_step_output=timedelta(minutes=60),
                   outfile='test_time_step60.nc')
        self.assertEqual(o1.history.shape, (100,25))
        self.assertEqual(o2.history.shape, (100,13))
        # Check that start and end conditions (longitudes) are idential
        self.assertItemsEqual(o1.history['lon'][:,24].compressed(), 
                              o2.history['lon'][:,12].compressed())
        self.assertItemsEqual(o1.history['lon'][:,0].compressed(),
                              o2.history['lon'][:,0].compressed())
        # Check that also run imported from file is identical
        o1i = OceanDrift(loglevel=20)
        o1i.io_import_file('test_time_step30.nc')
        o2i = OceanDrift(loglevel=20)
        o2i.io_import_file('test_time_step60.nc')
        os.remove('test_time_step30.nc')
        os.remove('test_time_step60.nc')
        self.assertItemsEqual(o2i.history['lon'][:,12].compressed(),
                              o2.history['lon'][:,12].compressed())
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

    def test_reader_boundary(self):
        # Check that the element outside reader coverage is
        # not deactivated if fallback value exist
        o = OceanDrift()
        nordic3d = reader_ROMS_native.Reader(o.test_data_folder() +
            '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        lon = [12.0, 12.0]
        lat = [70.0, 70.5]
        o.add_reader(nordic3d)
        o.fallback_values['land_binary_mask'] = 0
        o.seed_elements(lon, lat, number=2, radius=0,
                        time=nordic3d.start_time)
        o.run(steps=2, time_step=3600)
        self.assertEqual(o.num_elements_active(), 2)
        self.assertEqual(o.num_elements_deactivated(), 0)
        # Check that the outside element is deactivated,
        # if no fallback value exists
        o = OceanDrift()
        del o.fallback_values['x_sea_water_velocity']
        o.add_reader(nordic3d)
        o.fallback_values['land_binary_mask'] = 0
        o.seed_elements(lon, lat, number=2, radius=0,
                        time=nordic3d.start_time)
        o.run(steps=2, time_step=3600)
        self.assertEqual(o.num_elements_active(), 1)
        self.assertEqual(o.num_elements_deactivated(), 1)

    def test_reader_order(self):
        # Check that we get the same output indepenently of reader order
        o = OceanDrift(loglevel=50)
        norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
            '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
        basemap = reader_basemap_landmask.Reader(
            llcrnrlon=2, llcrnrlat=59.8, urcrnrlon=6, urcrnrlat=61,
            resolution='c', projection='merc')
        lon=4.; lat=60.

        # First run
        o.add_reader([basemap, norkyst, arome])
        o.seed_elements(lon, lat, time=norkyst.start_time)
        o.run(steps=30)
        # Second run
        # Check that we get almost identical results with other projection
        o1 = OceanDrift(loglevel=50)
        o1.add_reader([norkyst, arome, basemap])
        o1.seed_elements(lon, lat, time=norkyst.start_time)
        o1.run(steps=30)
        self.assertAlmostEqual(o.elements.lon, o1.elements.lon, 2)
        self.assertAlmostEqual(o.elements.lat, o1.elements.lat, 2)
        # Third run
        # Check that this is identical to run 1 if projection set equal
        o2 = OceanDrift(loglevel=50)
        o2.add_reader([norkyst, arome, basemap])
        o2.seed_elements(lon, lat, time=norkyst.start_time)
        o2.set_projection(basemap.proj4)
        o2.run(steps=30)
        self.assertEqual(o.elements.lon, o2.elements.lon)

    def test_seed_seafloor(self):
        o = OpenOil3D(loglevel=0)
        reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        # Landmask (Basemap)
        reader_basemap = reader_basemap_landmask.Reader(
                            llcrnrlon=4, llcrnrlat=61.0,
                            urcrnrlon=7, urcrnrlat=64,
                            resolution='c', projection='merc')
        o.add_reader([reader_basemap, reader_norkyst])
        lon = 4.5; lat = 62.0
        o.seed_elements(lon, lat, z='seafloor', time=reader_norkyst.start_time,
                        density=1000)
        o.config['processes']['turbulentmixing'] = True
        o.run(steps=3, time_step=900, time_step_output=900)
        #o.plot_property('z')
        z, status = o.get_property('z')
        self.assertAlmostEqual(z[0,0], -151.2, 1)  # Seeded at seafloor depth
        self.assertAlmostEqual(z[-1,0], -144, 2)  # After some rising

    def test_lift_above_seafloor(self):
        # See an element at some depth, and progapate towards coast
        # (shallower water) and check that it is not penetrating seafloor
        o = OceanDrift3D(loglevel=0, proj4='+proj=merc')
        reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
        reader_norkyst.buffer = 200
        o.add_reader([reader_norkyst],
                     variables='sea_floor_depth_below_sea_level')
        o.fallback_values['x_sea_water_velocity'] = 100 # Pure eastward motion
        o.fallback_values['y_sea_water_velocity'] = 0   
        o.fallback_values['land_binary_mask'] = 0   
        o.seed_elements(3.0, 62.0, z=-200, time=reader_norkyst.start_time)
        o.config['processes']['turbulentmixing'] = False
        o.run(steps=26, time_step=30)
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

if __name__ == '__main__':
    unittest.main()
