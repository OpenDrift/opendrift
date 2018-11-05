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


class TestPhysics(unittest.TestCase):
    """Tests for some physical parameterisations"""

    def test_droplet_diameters(self):
        o = OpenOil3D(loglevel=20, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 0
        o.fallback_values['y_wind'] = 0
        o.fallback_values['x_sea_water_velocity'] = 0
        o.fallback_values['y_sea_water_velocity'] = 0
        o.fallback_values['sea_surface_wave_period_at_variance_spectral_density_maximum'] = 5.8
        o.fallback_values['sea_surface_wave_significant_height'] = 3
        o.set_config('turbulentmixing:verticalresolution', 2)
        o.set_config('turbulentmixing:timestep', 4)
        # Setting droplet size range for subsea blowout
        o.set_config('input:spill:droplet_diameter_min_subsea', 0.0005)
        o.set_config('input:spill:droplet_diameter_max_subsea', 0.005)
        # Setting droplet size range for wave breaking
        o.set_config('turbulentmixing:droplet_diameter_min_wavebreaking', 1e-5)
        o.set_config('turbulentmixing:droplet_diameter_max_wavebreaking', 1e-3)
        # Number distribution frm Delvigne & Sweeney (s=-2.3):
        o.set_config('turbulentmixing:droplet_size_exponent', -2.3)
        o.seed_elements(4, 60, number=100, time=datetime.now(), z=-100)
        o.run(duration=timedelta(hours=3), time_step_output=900, time_step=900,
              stop_on_error=True)
        d_start = o.history['diameter'][:,0]
        d_end = o.history['diameter'][:,-1]
        # Check initial droplet sizes (expect range 0.0005 to 0.005)
        self.assertTrue(d_start.min() > o.get_config('input:spill:droplet_diameter_min_subsea'))
        self.assertTrue(d_start.max() < o.get_config('input:spill:droplet_diameter_max_subsea'))

    def test_constant_droplet_diameters(self):
        o = OpenOil3D(loglevel=50, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 0
        o.fallback_values['y_wind'] = 0
        o.fallback_values['x_sea_water_velocity'] = 0
        o.fallback_values['y_sea_water_velocity'] = 0
        o.fallback_values['sea_surface_wave_period_at_variance_spectral_density_maximum'] = 5.8
        o.fallback_values['sea_surface_wave_significant_height'] = 2.5
        o.set_config('turbulentmixing:verticalresolution', 2)
        o.set_config('turbulentmixing:timestep', 4)
        # Setting droplet size range for subsea blowout
        o.set_config('input:spill:droplet_diameter_min_subsea', 0.0005)
        o.set_config('input:spill:droplet_diameter_max_subsea', 0.005)
        o.set_config('turbulentmixing:droplet_size_exponent', -2.3)
        # Setting droplet size range for wave breaking
        o.set_config('turbulentmixing:droplet_diameter_min_wavebreaking', 1e-6)
        o.set_config('turbulentmixing:droplet_diameter_max_wavebreaking', 1e-3)
        diameter = 1e-4
        o.seed_elements(4, 60, number=100, time=datetime.now(),
                        diameter=diameter, z=-200)
        o.run(duration=timedelta(hours=2), time_step_output=900, time_step=900)

        d_start = o.history['diameter'][:,0]
        d_end = o.history['diameter'][:,-1]
        # Check droplet sizes before wavebreaking
        self.assertAlmostEqual(d_start.min(), diameter)
        self.assertAlmostEqual(d_start.max(), diameter)
        self.assertAlmostEqual(d_end.min(), diameter)
        self.assertAlmostEqual(d_end.min(), diameter)

    def test_vertical_mixing_nomixing(self):
        """Testing vertical mixing scheme

        The OpenOil3D model is run with various wind speeds,
        oil droplet diameters and time steps.
        The maximum depth of the particles is checked for sanity."""

        ########################################################
        # No wind/waves (i.e. no mixing)
        o = OpenOil3D(loglevel=20, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 0
        o.fallback_values['y_wind'] = 0
        o.fallback_values['x_sea_water_velocity'] = 0
        o.fallback_values['y_sea_water_velocity'] = 0
        o.seed_elements(4, 60, number=100, time=datetime.now())
        o.set_config('turbulentmixing:timestep', 5)
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
        #o.set_config('turbulentmixing:verticalresolution', 2)
        #o.set_config('turbulentmixing:timestep', 4)
        #o.run(duration=timedelta(hours=2), time_step_output=900, time_step=900)
        ##o.plot_property('z')
        ##o.plot_vertical_distribution()
        ##o.animation_profile()
        ## Check minimum depth
        #self.assertAlmostEqual(o.elements.z.min(), -46.0, 1)
        ########################################################

    def test_vertical_mixing_plantoil(self):
        #######################################################
        # 2.5m Hs, 10 mum radius (PlantOil)
        # Benchmark test from Jones et al. (2016)
        # NB: Entrainment length scale is not varied as in paper
        o = OpenOil3D(loglevel=30, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['sea_surface_wave_period_at_variance_spectral_density_maximum'] = 5.8
        o.fallback_values['sea_surface_wave_significant_height'] = 2.5
        o.fallback_values['x_wind'] = 10
        o.fallback_values['y_wind'] = 0
        o.fallback_values['x_sea_water_velocity'] = 0
        o.fallback_values['y_sea_water_velocity'] = 0
        o.seed_elements(4, 60, number=1000, diameter=0.00002,  # r = 10 micron
                        density=865, time=datetime.now(),
                        entrainment_length_scale=0.01)
        o.set_config('turbulentmixing:verticalresolution', 2)
        o.set_config('turbulentmixing:timestep', 4)
        o.run(duration=timedelta(hours=2), time_step_output=900, time_step=900)
        #o.plot_property('z')
        #o.plot_vertical_distribution()
        #o.animation_profile()
        # Check minimum depth
        self.assertAlmostEqual(o.elements.z.min(), -46.4, 1)
        #######################################################

    def test_vertical_mixing_plantoil_windonly(self):
        #######################################################
        # Same as above, but parameterising waves from wind
        o = OpenOil3D(loglevel=20, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 10
        o.fallback_values['y_wind'] = 0
        o.fallback_values['x_sea_water_velocity'] = 0
        o.fallback_values['y_sea_water_velocity'] = 0
        o.seed_elements(4, 60, number=1000, diameter=0.00002,  # r = 10 micron
                        density=865, time=datetime.now(),
                        entrainment_length_scale=0.01)
        o.set_config('turbulentmixing:verticalresolution', 2)
        o.set_config('turbulentmixing:timestep', 4)
        o.run(duration=timedelta(hours=2), time_step_output=900, time_step=900)
        #o.plot_vertical_distribution()
        self.assertAlmostEqual(o.elements.z.min(), -34.4, 1)
        #######################################################


        ########################################################
        ## Repeating last run, but with larger major (outer) time step
        ## Max mixing depth is expected to be same, but is slightly different
        ## This test is made to pass, but results should be checked
        o = OpenOil3D(loglevel=20, weathering_model='default')
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 10
        o.fallback_values['y_wind'] = 0
        o.fallback_values['x_sea_water_velocity'] = 0
        o.fallback_values['y_sea_water_velocity'] = 0
        o.seed_elements(4, 60, number=1000, diameter=0.00002,  # r = 10 micron
                        density=865, time=datetime.now(),
                        entrainment_length_scale=0.01)
        o.set_config('turbulentmixing:verticalresolution', 2)
        o.set_config('turbulentmixing:timestep', 4)
        o.run(duration=timedelta(hours=2),
              time_step_output=1800, time_step=1800)
        #o.plot_vertical_distribution()
        self.assertAlmostEqual(o.elements.z.min(), -34.1, 1)
        ########################################################


    def test_parameterised_stokes(self):
        o = OpenOil3D(loglevel=30, weathering_model='default')
        o.set_config('drift:use_tabularised_stokes_drift', False)
        o.set_config('processes:evaporation', False)
        o.fallback_values['land_binary_mask'] = 0
        o.fallback_values['x_wind'] = 10
        o.fallback_values['y_wind'] = 0
        o.fallback_values['x_sea_water_velocity'] = 0
        o.fallback_values['y_sea_water_velocity'] = 0
        o.seed_elements(lon=3, lat=60, time=datetime.now())
        o.run(steps=2)
        # Second run with parameterised Stokes drift
        o2 = OpenOil3D(loglevel=30, weathering_model='default')
        o2.set_config('drift:use_tabularised_stokes_drift', True)
        o2.set_config('processes:evaporation', False)
        o2.fallback_values['land_binary_mask'] = 0
        o2.fallback_values['x_wind'] = 10
        o2.fallback_values['y_wind'] = 0
        o2.fallback_values['x_sea_water_velocity'] = 0
        o2.fallback_values['y_sea_water_velocity'] = 0
        o2.seed_elements(lon=3, lat=60, time=datetime.now())
        o2.run(steps=2)
        # Check that stokes drift moves elements downwind
        self.assertTrue(o2.elements.lon > o.elements.lon)

if __name__ == '__main__':
    unittest.main()
