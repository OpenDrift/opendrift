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
        for oiltype in o.oiltypes[0:20]:
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

if __name__ == '__main__':
    unittest.main()
