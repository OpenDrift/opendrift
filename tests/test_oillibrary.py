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
import logging
from datetime import datetime, timedelta

import numpy as np

from opendrift.models.openoil3D import OpenOil3D

try:
    from oil_library import get_oil_props
    has_oil_library = True
except Exception as e:
    has_oil_library = False

class TestOil(unittest.TestCase):
    """Tests for OilLibrary"""

    @unittest.skipIf(has_oil_library is False,
                     'NOAA OilLibrary is needed')
    def test_oils(self):
        o = OpenOil3D(loglevel=50, weathering_model='noaa')
        for oiltype in o.oiltypes[10:15]:
            if oiltype == 'JP-8':
                continue
            o = OpenOil3D(loglevel=50, weathering_model='noaa')
            o.fallback_values['x_wind'] = 7
            o.fallback_values['x_sea_water_velocity'] = .7
            o.fallback_values['land_binary_mask'] = 0
            o.seed_elements(4.7, 60.0, radius=3000, number=3, z=0,
            time=datetime.now(), oiltype=oiltype)
            o.set_config('processes:evaporation', True)
            o.set_config('processes:emulsification', True)
            o.set_config('processes:turbulentmixing', False)
            o.run(steps=3)
            self.assertEqual(o.elements.mass_evaporated.min(),
                             o.elements.mass_evaporated.max())
            self.assertTrue(o.elements.mass_evaporated.min() > 0)
            self.assertTrue(o.elements.mass_evaporated.max() <= 1)
            print oiltype, o.elements.mass_evaporated.min()

    @unittest.skipIf(has_oil_library is False,
                     'NOAA OilLibrary is needed')
    def test_dispersion(self):
        for oil in ['SMORBUKK KONDENSAT', 'SKRUGARD']:
            for windspeed in [3, 8]:
                if oil == 'SKRUGARD' and windspeed == 3:
                    continue
                o = OpenOil3D(loglevel=50, weathering_model='noaa')
                o.seed_elements(lon=4.8, lat=60, number=100,
                                time=datetime.now(), oiltype=oil)
                o.set_config('processes:dispersion', True)
                o.set_config('wave_entrainment:droplet_size_distribution', 'Exponential')
                o.fallback_values['land_binary_mask'] = 0
                o.fallback_values['x_wind'] = windspeed
                o.fallback_values['y_sea_water_velocity'] = .3
                o.run(duration=timedelta(hours=6), time_step=900)

                b = o.get_oil_budget()
                actual_dispersed = b['mass_dispersed']/b['mass_total']
                print 'Dispersion fraction %f for ' \
                      '%s and wind speed %f' % \
                      (actual_dispersed[-1], oil, windspeed)
                if oil == 'SMORBUKK KONDENSAT' and windspeed == 3:
                    fraction_dispersed = 0
                elif oil == 'SMORBUKK KONDENSAT' and windspeed == 8:
                    fraction_dispersed = 0.11863789
                    fraction_dispersed = 0.1198852
                elif oil == 'SKRUGARD' and windspeed == 8:
                    fraction_dispersed =  0.23420876
                    fraction_dispersed =  0.23335621
                else:
                    fraction_dispersed = -1  # not defined
                self.assertAlmostEqual(actual_dispersed[-1],
                                       fraction_dispersed)
                #o.plot_oil_budget()

    @unittest.skipIf(has_oil_library is False,
                     'NOAA OilLibrary is needed')
    def test_no_dispersion(self):
        o = OpenOil3D(loglevel=50, weathering_model='noaa')
        o.seed_elements(lon=4.8, lat=60, number=100,
                        time=datetime.now(), oiltype='SKRUGARD')
        o.set_config('processes:dispersion', False)
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 8
        o.fallback_values['y_sea_water_velocity'] = .3
        o.run(duration=timedelta(hours=2), time_step=1800)
        b = o.get_oil_budget()
        actual_dispersed = b['mass_dispersed']/b['mass_total']
        self.assertAlmostEqual(actual_dispersed[-1], 0)

    @unittest.skipIf(has_oil_library is False,
                     'NOAA OilLibrary is needed')
    def test_droplet_distribution(self):
        for droplet_distribution in ['Johansen et al. (2015)',
                                     'Exponential']:
            o = OpenOil3D(loglevel=50, weathering_model='noaa')
            o.set_config('wave_entrainment:droplet_size_distribution',
                         droplet_distribution)
            o.seed_elements(lon=4.8, lat=60, number=100,
                            time=datetime.now(), oiltype='SKRUGARD')
            o.fallback_values['land_binary_mask'] = 0
            o.fallback_values['x_wind'] = 8
            o.fallback_values['y_sea_water_velocity'] = .3
            o.run(duration=timedelta(hours=1), time_step=1800)
            d = o.elements.diameter
            # Suspicious, Sintef-param should give larer droplets
            if droplet_distribution == 'Exponential':
                self.assertAlmostEqual(d.mean(), 0.0004933735)
            elif droplet_distribution == 'Johansen et al. (2015)':
                #self.assertAlmostEqual(d.mean(), 0.000072158)
                self.assertAlmostEqual(d.mean(), 0.00007332542)


if __name__ == '__main__':
    unittest.main()
