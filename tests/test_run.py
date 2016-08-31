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
from opendrift.models.openoil3D import OpenOil3D


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

    def test_output_time_step(self):
        o1 = OceanDrift(loglevel=20)
        norkyst = reader_netCDF_CF_generic.Reader(o1.test_data_folder() +
            '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        basemap = reader_basemap_landmask.Reader(
            llcrnrlon=4, llcrnrlat=59.8, urcrnrlon=6, urcrnrlat=61,
            resolution='h', projection='merc')
        o1.add_reader([basemap, norkyst])
        o1.seed_elements(5.25, 60.2, radius=3000, number=100,
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
        o2.seed_elements(5.25, 60.2, radius=3000, number=100,
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
                            resolution='h', projection='merc')
        o.add_reader([reader_basemap, reader_norkyst])
        lon = 4.5; lat = 62.0
        o.seed_elements(lon, lat, z='seafloor', time=reader_norkyst.start_time,
                        density=1000)
        o.config['processes']['turbulentmixing'] = True
        o.run(steps=3, time_step=900, time_step_output=900)
        #o.plot_property('z')
        z, status = o.get_property('z')
        self.assertAlmostEqual(z[0,0], -151.2, 1)  # Seeded at seafloor depth
        self.assertAlmostEqual(z[-1,0], -143, 2)  # After some rising


if __name__ == '__main__':
    unittest.main()
