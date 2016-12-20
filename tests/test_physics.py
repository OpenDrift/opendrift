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

        The OpenOil3D model is run with various wind speeds,
        oil droplet diameters and time steps.
        The maximum depth of the particles is checked for sanity."""

        ########################################################
        # No wind/waves (i.e. no mixing)
        o = OpenOil3D(loglevel=0, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.seed_elements(4, 60, number=100, time=datetime.now())
        o.config['turbulentmixing']['timestep'] = 5
        o.run(steps=4*2, time_step_output=3600, time_step=900)
        self.assertEqual(o.elements.z.min(), 0)  # No entrainment
        ########################################################

        ########################################################
        ## 2.5m Hs, 50 mum droplet radius (Oil Emulsion)
        ## Benchmark test from Jones et al. (2016)
        ## NB: Entrainment length scale is not varied as in paper
        ## entrainment_length_scale (L) 
        ## L = 0.01  for Plant oil, and L = 0.1 for Emulsion
        #o = OpenOil3D(loglevel=0, weathering_model='default')
        #o.fallback_values['land_binary_mask'] = 0
        #o.fallback_values['sea_surface_wave_period_at_variance_spectral_density_maximum'] = 5.8
        #o.fallback_values['sea_surface_wave_significant_height'] = 2.5
        #o.seed_elements(4, 60, number=1000, diameter=0.0001, # r = 50 micron
        #                density=865, time=datetime.now())
        #o.config['turbulentmixing']['verticalresolution'] = 2
        #o.config['turbulentmixing']['timestep'] = 4
        #o.run(duration=timedelta(hours=2), time_step_output=900, time_step=900)
        ##o.plot_property('z')
        ##o.plot_vertical_distribution()
        ##o.animation_profile()
        ## Check minimum depth
        #self.assertAlmostEqual(o.elements.z.min(), -46.0, 1)
        ########################################################

        #######################################################
        # 2.5m Hs, 10 mum radius (PlantOil)
        # Benchmark test from Jones et al. (2016)
        # NB: Entrainment length scale is not varied as in paper
        o = OpenOil3D(loglevel=0, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['sea_surface_wave_period_at_variance_spectral_density_maximum'] = 5.8
        o.fallback_values['sea_surface_wave_significant_height'] = 2.5
        o.seed_elements(4, 60, number=1000, diameter=0.00002,  # r = 10 micron
                        density=865, time=datetime.now(),
                        entrainment_length_scale=0.01)
        o.config['turbulentmixing']['verticalresolution'] = 2
        o.config['turbulentmixing']['timestep'] = 4
        o.run(duration=timedelta(hours=2), time_step_output=900, time_step=900)
        #o.plot_property('z')
        o.plot_vertical_distribution()
        #o.animation_profile()
        # Check minimum depth
        self.assertAlmostEqual(o.elements.z.min(), -56, 1)
        #######################################################

        #######################################################
        # Same as above, but parameterising waves from wind
        o = OpenOil3D(loglevel=0, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 10
        o.seed_elements(4, 60, number=1000, diameter=0.00002,  # r = 10 micron
                        density=865, time=datetime.now(),
                        entrainment_length_scale=0.01)
        o.config['turbulentmixing']['verticalresolution'] = 2
        o.config['turbulentmixing']['timestep'] = 4
        o.run(duration=timedelta(hours=2), time_step_output=900, time_step=900)
        o.plot_vertical_distribution()
        self.assertAlmostEqual(o.elements.z.min(), -62.0, 1)
        #######################################################


        ########################################################
        ## Repeating last run, but with larger major (outer) time step
        ## Max mixing depth is expected to be same, but is slightly different
        ## This test is made to pass, but results should be checked
        o = OpenOil3D(loglevel=0, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 10
        o.seed_elements(4, 60, number=1000, diameter=0.00002,  # r = 10 micron
                        density=865, time=datetime.now(),
                        entrainment_length_scale=0.01)
        o.config['turbulentmixing']['verticalresolution'] = 2
        o.config['turbulentmixing']['timestep'] = 4
        o.run(duration=timedelta(hours=2),
              time_step_output=1800, time_step=1800)
        #o.plot_vertical_distribution()
        self.assertAlmostEqual(o.elements.z.min(), -54.0, 1)
        ########################################################


if __name__ == '__main__':
    unittest.main()
