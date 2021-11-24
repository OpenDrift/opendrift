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
import numpy as np

from opendrift.models.oceandrift import OceanDrift
from opendrift.models.openoil import OpenOil

class TestSeed(unittest.TestCase):
    """Tests of seeding methods"""

    def test_seed_elements(self):
        # Some cases with expected outcome
        lon_vec = np.array([3, 4, 5])
        lat_vec = np.array([63, 64, 65])
        lon0 = 3
        lat0 = 60
        t0 = datetime.now()
        t1 = t0 + timedelta(hours=6)
        t2 = t1 + timedelta(hours=6)
        t_vec = [t0, t1, t2]
        number_config = 222

        cases = [
            {'lon': lon0, 'lat': lat0, 'time': t0,
                  'number': None, 'expected': number_config},
            {'lon': lon0, 'lat': lat0, 'time': t0,
                  'number': 12, 'expected': 12, 'maxlat': lat0},
            {'lon': lon0, 'lat': lat0, 'time': t0, 'radius': 1000,
                  'number': 12, 'expected': 12, 'maxlat': lat0 + .013},
            {'lon': lon0, 'lat': lat0, 'time': [t0, t1],
                  'maxtime': t1, 'number': 12, 'expected': 12},
            {'lon': lon0, 'lat': lat0, 'time': t_vec, # Time series
                  'maxtime': t2, 'number': None, 'expected': len(t_vec)},
            {'lon': lon0, 'lat': lat0, 'time': t_vec,
                  'number': 12, 'expected': 'error'},
            {'lon': lon_vec, 'lat': lat_vec, 'time': t0,
                  'number': None, 'expected': 3},
            {'lon': lon_vec, 'lat': lat_vec, 'time': t0,
                  'wind_drift_factor': .05,
                  'number': None, 'expected': 3},
            {'lon': lon_vec, 'lat': lat_vec, 'time': t0,
                  'wind_drift_factor': [.01, .02, .03],
                  'number': None, 'expected': 3},
            {'lon': lon_vec, 'lat': lat_vec, 'time': t0,
                  'number': 4, 'expected': 'error'},
            {'lon': lon_vec, 'lat': lat_vec, 'time': t_vec,
                  'number': None, 'expected': 3, 'maxtime': t_vec[-1]},
        ]

        for case in cases:
            o = OceanDrift(loglevel=50)
            o._set_config_default('seed:number', number_config)
            expected = case['expected']
            del case['expected']
            if 'maxlat' in case:
                maxlat = case['maxlat']
                del case['maxlat']
            else:
                maxlat = None
            if 'maxtime' in case:
                maxtime = case['maxtime']
                del case['maxtime']
            else:
                maxtime = None
            if expected == 'error':
                with self.assertRaises(ValueError):
                    n = o.seed_elements(**case)
            else:
                o.seed_elements(**case)
                self.assertEqual(o.num_elements_total(), expected)
                if maxlat is not None:
                    self.assertAlmostEqual(o.elements_scheduled.lat.max(),
                                           maxlat, 2)
                if maxtime is not None:
                    self.assertEqual(
                        o.elements_scheduled_time.max().replace(
                            microsecond=0),
                        maxtime.replace(microsecond=0))
                if 'wind_drift_factor' in case:
                    self.assertEqual(
                        np.array(case['wind_drift_factor']
                                ).max().astype(np.float32),
                        o.elements_scheduled.wind_drift_factor.max())

    def test_seed_cone(self):

        # Some cases with expected outcome
        lon0 = 3
        lon1 = 4
        lon2 = 5
        lat0 = 60
        lat1 = 60.2
        lat2 = 60.4
        lon_vec = np.array([lon0, lon1, lon2])
        lat_vec = np.array([lat0, lat1, lat2])
        t0 = datetime.now()
        t1 = t0 + timedelta(hours=6)
        t2 = t1 + timedelta(hours=6)
        t_vec = [t0, t1, t2]
        r0 = 1000
        r1 = 2000
        number_config = 222

        cases = [
            {'lon': lon0, 'lat': lat0, 'time': t0,
                  'number': None, 'expected': number_config},
            {'lon': lon0, 'lat': lat0, 'time': t0,
                  'number': 12, 'expected': 12, 'maxlat': lat0},
            {'lon': lon0, 'lat': lat0, 'time': t0, 'radius': 1000,
                  'number': 12, 'expected': 12, 'maxlat': lat0 + .013},
            {'lon': lon0, 'lat': lat0, 'time': [t0, t1],
                  'maxtime': t1, 'number': 12, 'expected': 12},
            {'lon': lon0, 'lat': lat0, 'time': t_vec,
                  'number': 12, 'expected': 'error'},
            {'lon': lon_vec, 'lat': lat_vec, 'time': t0,
                  'number': None, 'expected': 'error'},
            {'lon': lon0, 'lat': lat0, 'time': t0,
                  'wind_drift_factor': .05,
                  'number': None, 'expected': number_config},
            {'lon': [lon0, lon1], 'lat': [lat0, lat1], 'time': t0,
                  'number': None, 'expected': number_config},
            {'lon': [lon0, lon1], 'lat': [lat0, lat1], 'time': t0,
                  'number': 12, 'expected': 12},
            {'lon': [lon0, lon1], 'lat': [lat0, lat1], 'time': t0,
                    'radius': 0, 'number': 200, 'expected': 200, 'maxlat': lat1},
            {'lon': [lon0, lon1], 'lat': [lat0, lat1], 'time': t0,
                    'radius': 2000, 'number': 200, 'expected': 200, 'maxlat': lat1 + .024},
            {'lon': [lon0, lon1], 'lat': [lat0, lat1], 'time': t0,
                    'radius': [1000, 2000], 'number': 200, 'expected': 200, 'maxlat': lat1 + .024},
            {'lon': [lon0, lon1], 'lat': [lat0, lat1], 'time': t0,
                    'radius': [1000, 2000, 3000], 'number': 200, 'expected': 'error'},
                  #'wind_drift_factor': [.01, .02, .03],
            {'lon': [lon0, lon1], 'lat': [lat0, lat1], 'time': [t0, t1],
                  'number': None, 'expected': number_config, 'maxtime': t1},
            {'lon': [lon0, lon1], 'lat': [lat0, lat1], 'time': [t0, t1], 'radius': [r0, r1],
                  'number': None, 'expected': number_config, 'maxtime': t1},
            {'lon': [lon0, lon1], 'lat': [lat0, lat1], 'time': [t0, t1], 'radius': [r0, r0],
                  'number': None, 'expected': number_config, 'maxtime': t1},
        ]

        for case in cases:
            o = OceanDrift(loglevel=50)
            o._set_config_default('seed:number', number_config)
            expected = case['expected']
            del case['expected']
            if 'maxlat' in case:
                maxlat = case['maxlat']
                del case['maxlat']
            else:
                maxlat = None
            if 'maxtime' in case:
                maxtime = case['maxtime']
                del case['maxtime']
            else:
                maxtime = None
            if expected == 'error':
                with self.assertRaises(ValueError):
                    n = o.seed_cone(**case)
            else:
                o.seed_cone(**case)
                self.assertEqual(o.num_elements_total(), expected)
                if maxlat is not None:
                    self.assertAlmostEqual(o.elements_scheduled.lat.max(),
                                           maxlat, 2)
                if maxtime is not None:
                    self.assertEqual(
                        o.elements_scheduled_time.max().replace(
                            microsecond=0),
                        maxtime.replace(microsecond=0))
                if 'wind_drift_factor' in case:
                    self.assertEqual(
                        np.array(case['wind_drift_factor']
                                ).max().astype(np.float32),
                        o.elements_scheduled.wind_drift_factor.max())

    def test_seed_seafloor(self):
        # Check that seed:seafloor overrides z to 'seafloor'
        o = OceanDrift(loglevel=50)
        from opendrift.readers import reader_constant
        r = reader_constant.Reader({
            'sea_floor_depth_below_sea_level': 200})
        o.add_reader(r)
        o.set_config('seed:seafloor', True)
        o.seed_elements(lon=4, lat=60, time=datetime.now())
        self.assertAlmostEqual(o.elements_scheduled.z[0], -200)

    def test_seed_ocean_only(self):
        o = OceanDrift(loglevel=50)
        o.set_config('seed:ocean_only', True)
        o.seed_elements(lon=[4.8, 5.15], lat=[60, 60], time=datetime.now())
        o.run(steps=1)
        lat, s = o.get_property('lat')[0]
        # Check that the second point over land has been moved to ocean
        self.assertAlmostEqual(lat[0], 60, 5)
        self.assertNotAlmostEqual(lat[1], 60, 5)

    def test_seed_letters(self):
        o = OceanDrift(loglevel=50)
        o.seed_letters('Obey Soros', lon=-2, lat=61, time=datetime.now(), number=1000)
        self.assertAlmostEqual(o.elements_scheduled.lon.max(), 5.64, 1)

if __name__ == '__main__':
    unittest.main()
