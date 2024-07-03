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
import numpy as np

from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_oscillating
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.pelagicegg import PelagicEggDrift
from opendrift.models.oceandrift import OceanDrift


class TestStranding(unittest.TestCase):

    def test_stranding_3d(self):
        o = PelagicEggDrift(loglevel=30)
        reader_nordic = reader_ROMS_native.Reader(o.test_data_folder() +
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        o.add_reader(reader_nordic)
        o.set_config('environment:fallback:y_wind', 10)  # Some wind for mixing
        o.set_config('general:coastline_action', 'stranding')
        o.set_config('general:coastline_approximation_precision', None)
        o.set_config('vertical_mixing:timestep', 120)
        o.set_config('drift:max_speed', .1)
        o.seed_elements(lon=14.0, lat=68.15, radius=2000, number=100,
                        time=[reader_nordic.start_time,
                              reader_nordic.end_time], z=0)

        o.run(end_time=reader_nordic.end_time, time_step=3600*6)
        self.assertEqual(o.status_categories[1], 'stranded')
        self.assertEqual(o.elements_deactivated.status.min(), 1)
        self.assertEqual(o.elements_deactivated.status.max(), 1)
        self.assertEqual(o.num_elements_scheduled(), 0)
        self.assertEqual(o.num_elements_active(), 82)
        self.assertEqual(o.num_elements_activated(), 100)
        self.assertEqual(o.num_elements_deactivated(), 18)
        self.assertEqual(o.num_elements_total(), 100)

        # Check calculated trajectory lengths and speeds
        total_length, distances, speeds = o.get_trajectory_lengths()
        self.assertAlmostEqual(total_length.max(), 15722.7, 1)
        self.assertAlmostEqual(total_length.min(), 1225.2, 1)
        self.assertAlmostEqual(speeds.max(), 0.132, 1)
        self.assertAlmostEqual(distances.max(), 2859.0, 1)

    def test_stranding_roms(self):
        o = PelagicEggDrift(loglevel=20)
        o.set_config('environment:fallback:x_sea_water_velocity', 1)
        o.set_config('general:coastline_action', 'previous')
        o.set_config('drift:vertical_mixing', False)
        o.set_config('drift:max_speed', 1)
        reader_arctic = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
        '2Feb2016_Nordic_sigma_3d/Arctic20_1to5Feb_2016.nc')
        reader_nordic = reader_ROMS_native.Reader(o.test_data_folder() +
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        o.add_reader(reader_arctic)
        o.add_reader(reader_nordic)
        o.seed_elements(lon=13.0, lat=68.0, radius=20000, number=100,
                        time=[reader_arctic.start_time,
                              reader_nordic.end_time], z=-30)
        o.run(end_time=reader_nordic.end_time, time_step=3600*36)
        self.assertEqual(o.num_elements_scheduled(), 0)
        self.assertEqual(o.num_elements_active(), 100)
        self.assertEqual(o.num_elements_activated(), 100)
        self.assertEqual(o.num_elements_deactivated(), 0)
        self.assertEqual(o.num_elements_total(), 100)

    @pytest.mark.slow
    def test_stranding_options(self):
        reader_osc = reader_oscillating.Reader(
                'x_sea_water_velocity', amplitude=.5,
                period_seconds=3600*6,
                zero_time=datetime.now())

        reader_global = reader_global_landmask.Reader()

        # Three different stranding options, with
        # expected final status and position
        options = ['stranding', 'previous', 'none']
        status = ['stranded', 'active', 'active']
        lons_backward = [16.157, 16.087, 16.167]
        lons_forward = [16.198, 16.115, 16.167]

        for i, option in enumerate(options):
            for direction in ['forward', 'backward']:
                if direction == 'forward':
                    lons = lons_forward
                    time_step = 900
                else:
                    lons = lons_backward
                    time_step = -900
                o = OceanDrift(loglevel=50)
                o.set_config('general:coastline_action', option)
                o.set_config('general:coastline_approximation_precision', None)
                o.add_reader([reader_osc, reader_global])
                # Adding northwards drift
                o.set_config('environment:constant:y_sea_water_velocity', .1)
                #o.seed_elements(lon=13.35, lat=68.0, radius=0,
                o.seed_elements(lon=16.12, lat=68.5, radius=0,
                                time=reader_osc.zero_time)
                o.run(duration=timedelta(hours=22), time_step=time_step)
                #o.animation()
                #o.plot ()
                print('Testing stranding: %s, %s' % (option, direction))
                if len(o.elements) == 1:
                    el = o.elements
                else:
                    el = o.elements_deactivated
                self.assertEqual(o.status_categories[int(el.status)], status[i])
                self.assertIsNone(np.testing.assert_array_almost_equal(
                    el.lon, lons[i], 3))

    def test_interact_coastline_global(self):
        reader_global = reader_global_landmask.Reader()

        o = OceanDrift(loglevel=20)
        o.add_reader(reader_global)
        o.set_config('general:coastline_action', 'previous')
        o.set_config('general:use_auto_landmask', False)
        o.set_config('environment:fallback:x_sea_water_velocity', .7)
        o.seed_elements(lon=5, lat=60.49, time=datetime.now())
        o.run(time_step=3600, steps=30)
        lons = o.history['lon'][0]
        self.assertAlmostEqual(lons[0], 5, 2)
        self.assertAlmostEqual(lons[-2], 5.092, 2)
        self.assertAlmostEqual(lons[-1], 5.092, 2)

    def test_stranding_approximation(self):
        o = OceanDrift(loglevel=0)
        o.set_config('environment:constant:x_sea_water_velocity', 1)
        o.set_config('general:coastline_approximation_precision', None)
        o.seed_elements(lon=4.55, lat=60, time=datetime.now())
        o.run(steps=10)
        o2 = OceanDrift(loglevel=0)
        o2.set_config('environment:constant:x_sea_water_velocity', 1)
        o2.set_config('general:coastline_approximation_precision', .001)
        o2.seed_elements(lon=4.55, lat=60, time=datetime.now())
        o2.run(steps=10)
        self.assertAlmostEqual(o.elements_deactivated.lon[0], 5.066, 3)
        self.assertAlmostEqual(o2.elements_deactivated.lon[0], 5.051, 3)


if __name__ == '__main__':
    unittest.main()
