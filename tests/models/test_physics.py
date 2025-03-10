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
from datetime import datetime, timedelta
import os
import glob

import numpy as np
import trajan

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.models.openoil import OpenOil
from opendrift.models.physics_methods import verticaldiffusivity_Large1994, verticaldiffusivity_Sundby1983


class TestPhysics(unittest.TestCase):
    """Tests for some physical parameterisations"""

    def test_vertical_diffusivity(self):
        windspeeds = np.arange(0, 20, 5)
        depths = np.arange(0, 80, 5)
        wind, depth = np.meshgrid(windspeeds, depths)
        KLarge = verticaldiffusivity_Large1994(wind, depth)
        KSundby = verticaldiffusivity_Sundby1983(wind, depth)

        self.assertAlmostEqual(KLarge.min(), 0, 3)
        self.assertAlmostEqual(KLarge.max(), 0.2017, 3)
        self.assertAlmostEqual(KSundby.min(), 0, 3)
        self.assertAlmostEqual(KSundby.max(), 0.0585, 3)

    def test_droplet_diameters(self):
        o = OpenOil(loglevel=20)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.set_config('environment:fallback:sea_surface_wave_period_at_variance_spectral_density_maximum', 5.8)
        o.set_config('environment:fallback:sea_surface_wave_significant_height', 3)

        o.set_config('vertical_mixing:timestep', 4)
        # Setting droplet size range for subsea blowout
        o.set_config('seed:droplet_diameter_min_subsea', 0.0005)
        o.set_config('seed:droplet_diameter_max_subsea', 0.005)
        # Setting droplet size range for wave breaking
        o.seed_elements(4, 60, number=100, time=datetime.now(), z=-100)
        o.run(duration=timedelta(hours=3), time_step=900)
        d_start = o.result.diameter[:,0].values
        d_end = o.result.diameter[:,-1].values
        # Check initial droplet sizes (expect range 0.0005 to 0.005)
        self.assertTrue(d_start.min() >
                o.get_config('seed:droplet_diameter_min_subsea'))
        self.assertTrue(d_start.max() <
                o.get_config('seed:droplet_diameter_max_subsea'))

    def test_constant_droplet_diameters(self):
        o = OpenOil(loglevel=50)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.set_config('environment:fallback:sea_surface_wave_period_at_variance_spectral_density_maximum', 5.8)
        o.set_config('environment:fallback:sea_surface_wave_significant_height', 2.5)
        o.set_config('vertical_mixing:timestep', 4)
        # Setting droplet size range for subsea blowout
        o.set_config('seed:droplet_diameter_min_subsea', 0.0005)
        o.set_config('seed:droplet_diameter_max_subsea', 0.005)
        diameter = 1e-4
        o.seed_elements(4, 60, number=100, time=datetime.now(),
                        diameter=diameter, z=-200)
        o.run(duration=timedelta(hours=2), time_step_output=900, time_step=900)

        d_start = o.result.diameter[:,0].values
        d_end = o.result.diameter[:,-1].values
        # Check droplet sizes before wavebreaking
        self.assertAlmostEqual(d_start.min(), diameter)
        self.assertAlmostEqual(d_start.max(), diameter)
        self.assertAlmostEqual(d_end.min(), diameter)
        self.assertAlmostEqual(d_end.min(), diameter)

    def test_vertical_mixing_nomixing(self):
        """Testing vertical mixing scheme

        The OpenOil model is run with various wind speeds,
        oil droplet diameters and time steps.
        The maximum depth of the particles is checked for sanity."""

        ########################################################
        # No wind/waves (i.e. no mixing)
        o = OpenOil(loglevel=20)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 0)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.set_config('vertical_mixing:timestep', 5)
        o.seed_elements(4, 60, number=100, time=datetime.now())
        o.run(steps=4*2, time_step_output=3600, time_step=900)
        self.assertEqual(o.elements.z.min(), 0)  # No entrainment
        ########################################################

    def test_vertical_mixing_plantoil(self):
        #######################################################
        # 2.5m Hs, 10 mum radius (PlantOil)
        # Benchmark test from Jones et al. (2016)
        # NB: Entrainment length scale is not varied as in paper
        o = OpenOil(loglevel=30)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:sea_surface_wave_period_at_variance_spectral_density_maximum', 5.8)
        o.set_config('environment:fallback:sea_surface_wave_significant_height', 2.5)
        o.set_config('environment:fallback:x_wind', 10)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.set_config('vertical_mixing:timestep', 4)
        o.seed_elements(4, 60, number=1000, diameter=0.00002,  # r = 10 micron
                        density=865, time=datetime.now(), oil_type='AASGARD A 2003')
        o.run(duration=timedelta(hours=2), time_step_output=900, time_step=900)
        #o.plot_property('z')
        #o.plot_vertical_distribution()
        #o.animation_profile()
        # Check minimum depth
        self.assertAlmostEqual(o.elements.z.min(), -49.65, 1)
        #######################################################

    def test_vertical_mixing_plantoil_windonly(self):
        #######################################################
        # Same as above, but parameterising waves from wind
        o = OpenOil(loglevel=20)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 10)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.set_config('vertical_mixing:timestep', 4)
        o.seed_elements(4, 60, number=1000, diameter=0.00002,  # r = 10 micron
                        density=865, time=datetime.now(), oil_type='AASGARD A 2003')

        o.run(duration=timedelta(hours=2), time_step_output=900, time_step=900)
        #o.plot_vertical_distribution()
        self.assertAlmostEqual(o.elements.z.min(), -48.61, 1)
        #######################################################


        ########################################################
        ## Repeating last run, but with larger major (outer) time step
        ## Max mixing depth is expected to be same, but is slightly different
        ## This test is made to pass, but results should be checked
        o = OpenOil(loglevel=20)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 10)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.set_config('vertical_mixing:timestep', 4)
        o.seed_elements(4, 60, number=1000, diameter=0.00002,  # r = 10 micron
                        density=865, time=datetime.now(), oil_type='AASGARD A 2003')

        o.run(duration=timedelta(hours=2),
              time_step_output=1800, time_step=1800)
        #o.plot_vertical_distribution()
        self.assertAlmostEqual(o.elements.z.min(), -49.2, 1)
        ########################################################

    def test_verticalmixing_schemes(self):

        for scheme in ['environment', 'windspeed_Large1994',
                       'windspeed_Sundby1983', 'constant']:
            o = OpenOil(loglevel=50, weathering_model='noaa')
            o.set_config('environment:fallback:land_binary_mask', 0)
            o.set_config('environment:fallback:x_wind', 10)
            o.set_config('environment:fallback:y_wind', 0)
            o.set_config('environment:fallback:x_sea_water_velocity', 0)
            o.set_config('environment:fallback:y_sea_water_velocity', 0)
            o.set_config('environment:fallback:ocean_vertical_diffusivity', 0)
            #o.set_config('vertical_mixing:timestep', 4)
            o.set_config('vertical_mixing:diffusivitymodel', scheme)
            o.seed_elements(4, 60, number=1000, diameter=0.00002,  # r = 10 micron
                            density=865, time=datetime.now(), oil_type='AASGARD A 2003')

            o.run(duration=timedelta(hours=2), time_step=900)

            if scheme == 'environment':  # presently this is fallback
                self.assertAlmostEqual(o.elements.z.min(), -48.9, 1)
            elif scheme == 'windspeed_Large1994':
                self.assertAlmostEqual(o.elements.z.min(), -48.9, 1)
            elif scheme == 'windspeed_Sundby1983':
                self.assertAlmostEqual(o.elements.z.min(), -51.75, 1)
            elif scheme == 'constant':
                self.assertAlmostEqual(o.elements.z.min(), -3.57, 1)

    def test_parameterised_stokes(self):
        o = OpenOil(loglevel=30)
        o.set_config('drift:use_tabularised_stokes_drift', False)
        o.set_config('processes:evaporation', False)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 10)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.seed_elements(lon=3, lat=60, time=datetime.now())
        o.run(steps=2)
        # Second run with parameterised Stokes drift
        o2 = OpenOil(loglevel=30)
        o2.set_config('drift:use_tabularised_stokes_drift', True)
        o2.set_config('processes:evaporation', False)
        o2.set_config('environment:fallback:land_binary_mask', 0)
        o2.set_config('environment:fallback:x_wind', 10)
        o2.set_config('environment:fallback:y_wind', 0)
        o2.set_config('environment:fallback:x_sea_water_velocity', 0)
        o2.set_config('environment:fallback:y_sea_water_velocity', 0)
        o2.seed_elements(lon=3, lat=60, time=datetime.now())
        o2.run(steps=2)
        # Check that stokes drift moves elements downwind
        self.assertTrue(o2.elements.lon > o.elements.lon)
        
        
if __name__ == '__main__':
    unittest.main()
