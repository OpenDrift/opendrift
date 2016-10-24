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
import glob

import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.models.openoil3D import OpenOil3D


class TestRun(unittest.TestCase):
    """Tests for some physical parameterisations"""

    def test_vertical_mixing(self):
        """Testing vertical mixing scheme

        The OpenOil3D model is run for one hour with various
        wind speeds, oil droplet diameters and time steps.
        The maximum depth of the particles is checked for sanity."""

        #######################################################
        # No wind/waves (i.e. no mixing)
        o = OpenOil3D(loglevel=0, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.seed_elements(4, 60, number=100, time=datetime.now())
        o.config['turbulentmixing']['timestep'] = 5
        o.run(steps=4*2, time_step_output=3600, time_step=900)
        self.assertEqual(o.elements.z.min(), 0)  # No entrainment
        #######################################################

        #######################################################
        # 10 m/s wind speed, 0.001m droplet radius
        o = OpenOil3D(loglevel=0, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 10
        o.seed_elements(4, 60, number=100, diameter=.0002, # 200 micron
                        time=datetime.now())
        o.config['turbulentmixing']['timestep'] = 5
        o.run(steps=4*2, time_step_output=3600, time_step=900)
        # Droplets mixed down to 8 m depth
        self.assertAlmostEqual(o.elements.z.min(), -8.1, 1)
        #######################################################

        #######################################################
        # 10 m/s wind speed, 0.00001m droplet radius
        o = OpenOil3D(loglevel=0, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 10
        o.seed_elements(4, 60, number=100, diameter=.000002,  # 20 micron
                        time=datetime.now())
        o.config['turbulentmixing']['timestep'] = 5
        o.run(steps=4*2, time_step_output=3600, time_step=900)
        # Droplets mixed down to 42 m depth
        self.assertAlmostEqual(o.elements.z.min(), -42.0, 1)
        #######################################################

        #######################################################
        # Repeating last run, but with larger mixing time step (20 s vs 5s)
        # Max mixing depth is quite different, though no warning is issued
        # This test is made to pass, but results should be checked
        o = OpenOil3D(loglevel=0, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 10
        o.seed_elements(4, 60, number=100, diameter=.000002,  # 20 micron
                        time=datetime.now())
        o.config['turbulentmixing']['timestep'] = 20
        o.run(steps=4*2, time_step_output=3600, time_step=900)
        # Droplets mixed down to 51 m depth, due to larger mixing time step
        # This is not expected to make a difference, and should be checked!
        self.assertAlmostEqual(o.elements.z.min(), -51.0, 1)
        #######################################################

        #######################################################
        # Repeating last run, but with larger major (outer) time step
        o = OpenOil3D(loglevel=0, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 10
        o.seed_elements(4, 60, number=100, diameter=.000002,  # 20 micron
                        time=datetime.now())
        o.config['turbulentmixing']['timestep'] = 20
        o.run(steps=4*2, time_step_output=3600, time_step=3600)
        # Droplets mixed down to 75 m depth, due to larger major time step
        # This is not expected to make a difference, and should be checked!
        self.assertAlmostEqual(o.elements.z.min(), -75.0, 1)
        #######################################################


if __name__ == '__main__':
    unittest.main()
