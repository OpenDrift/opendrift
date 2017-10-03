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

from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_oscillating
from opendrift.models.pelagicegg import PelagicEggDrift
from opendrift.models.oceandrift import OceanDrift


class TestStranding(unittest.TestCase):

    def test_stranding_3d(self):
        o = PelagicEggDrift(loglevel=30)
        reader_nordic = reader_ROMS_native.Reader(o.test_data_folder() +
        '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
        o.add_reader(reader_nordic)
        o.fallback_values['y_wind'] = 10  # Some wind for mixing
        o.seed_elements(lon=14.0, lat=68.15, radius=2000, number=100,
                        time=[reader_nordic.start_time,
                              reader_nordic.end_time], z=0)
        o.set_config('general:basemap_resolution', 'i')
        o.set_config('general:coastline_action', 'stranding')
        o.set_config('turbulentmixing:timestep', 120)
        o.set_config('turbulentmixing:verticalresolution', 10)
        o.max_speed=.1
        o.run(end_time=reader_nordic.end_time, time_step=3600*6)
        self.assertEqual(o.status_categories[1], 'stranded')
        self.assertEqual(o.elements_deactivated.status.min(), 1)
        self.assertEqual(o.elements_deactivated.status.max(), 1)
        self.assertEqual(o.num_elements_scheduled(), 0)
        self.assertEqual(o.num_elements_active(), 86)
        self.assertEqual(o.num_elements_activated(), 100)
        self.assertEqual(o.num_elements_deactivated(), 14)
        self.assertEqual(o.num_elements_total(), 100)

    def test_stranding_options(self):

        reader_osc = reader_oscillating.Reader(
                'x_sea_water_velocity', amplitude=1,
                zero_time=datetime.now())

        reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=12, llcrnrlat=67.6,
                    urcrnrlon=13.6, urcrnrlat=68.1,
                    resolution='i', projection='merc')
 
        # Three different stranding options, with
        # expected final status and position
        options = ['stranding', 'previous', 'none']
        status = ['stranded', 'active', 'active']
        lons = [12.930, 13.348, 12.444]

        for i, option in enumerate(options):
            o = OceanDrift(loglevel=30)
            o.set_config('general:coastline_action', option)
            o.add_reader([reader_osc, reader_basemap])
            # Adding northwards drift
            o.fallback_values['y_sea_water_velocity'] = .2
            o.seed_elements(lon=12.2, lat=67.7, radius=0,
                            time=reader_osc.zero_time)
            o.run(steps=28, time_step=3600*2)
            print 'Testing stranding: %s' % option
            if len(o.elements) == 1:
                el = o.elements
            else:
                el = o.elements_deactivated
            self.assertEqual(o.status_categories[int(el.status)], status[i])
            self.assertAlmostEqual(el.lon, lons[i], 2)
            #o.plot()


if __name__ == '__main__':
    unittest.main()
