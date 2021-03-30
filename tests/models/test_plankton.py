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
import numpy as np
import pyproj

from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_ArtificialOceanEddy
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_constant
from opendrift.models.pelagicplankton_moana import PelagicPlanktonDrift


import opendrift
print(opendrift.versions())

class TestModels(unittest.TestCase):
    """Tests for OpenDrift models"""

    def test_PelagicPlankton_vertical_motion_pelagicegg(self):

        o = PelagicPlanktonDrift(loglevel=0)
        #s.list_configspec()
        o.set_config('environment:fallback:land_binary_mask', 0)

        # default values for forcing variables
        o.set_config('environment:fallback:x_sea_water_velocity', 1.0) #m/s 
        o.set_config('environment:fallback:y_sea_water_velocity', 0) #m/s 
        o.set_config('environment:fallback:x_wind', 0 )#m/s
        o.set_config('environment:fallback:y_wind', 0) #m/s
        o.set_config('environment:fallback:sea_water_temperature', 12.0 )# degrees Celsius
        o.set_config('environment:fallback:sea_water_salinity', 35.0 )   # ppt
        o.set_config('environment:fallback:ocean_vertical_diffusivity', 0.000) # m2/s

        # use the update_terminal_velocity_pelagicegg() for the vertical motion
        o.set_config('biology:vertical_migration_speed_constant',None) #'float(min=0.0, max=1e-3, default=None)', comment=' Constant vertical migration rate (m/s), if None, use values from update_terminal_velocity_pelagicegg()') #


        o.seed_elements(lon=2, lat=60,z = -1.0, time=datetime(2021, 2, 24, 3, 46, 8), number=1)
        o.run(duration=timedelta(hours=4))
        #self.assertAlmostEqual(s.elements.lon.max(), 2.1273, 3)  # Without setting config
        self.assertAlmostEqual(o.elements.lon.max(), 2.258064243, 3)     
        self.assertAlmostEqual(o.elements.z, -0.01, 3)  

    def test_PelagicPlankton_vertical_motion_constant(self):

        o = PelagicPlanktonDrift(loglevel=0)
        #s.list_configspec()
        o.set_config('environment:fallback:land_binary_mask', 0)

        # default values for forcing variables
        o.set_config('environment:fallback:x_sea_water_velocity', 1.0) #m/s 
        o.set_config('environment:fallback:y_sea_water_velocity', 0) #m/s 
        o.set_config('environment:fallback:x_wind', 0 )#m/s
        o.set_config('environment:fallback:y_wind', 0) #m/s
        o.set_config('environment:fallback:sea_water_temperature', 12.0 )# degrees Celsius
        o.set_config('environment:fallback:sea_water_salinity', 35.0 )   # ppt
        o.set_config('environment:fallback:ocean_vertical_diffusivity', 0.000) # m2/s

        # use the update_terminal_velocity_constant() for the vertical motion, with daytime and nighttime positions
        o.set_config('biology:vertical_migration_speed_constant',1.0e-4) #'float(min=0.0, max=1e-3, default=None)', comment=' Constant vertical migration rate (m/s), if None, use values from update_terminal_velocity_pelagicegg()') #
        o.set_config('biology:vertical_position_daytime', -5.0)#'float(min=-1000.0, max=0.0, default=-5.0)',   comment='the depth a species is expected to inhabit during the day time, in meters, negative down') #
        o.set_config('biology:vertical_position_nighttime', -1.0) #'float(min=-1000.0, max=0.0, default=-1.0)', comment='the depth a species is expected to inhabit during the night time, in meters, negative down') #

        o.seed_elements(lon=2, lat=60,z = -1.0, time=datetime(2021, 2, 24, 3, 46, 8), number=1)
        o.run(duration=timedelta(hours=4))
        
        self.assertAlmostEqual(o.elements.lon.max(), 2.258064243, 3)    
        self.assertAlmostEqual(o.elements.z.max(), -2.44, 3)  # note: different vertical position at end of the simulation that previous test


    def test_PelagicPlankton_survival_rate(self):

        o = PelagicPlanktonDrift(loglevel=0)
        #s.list_configspec()
        o.set_config('environment:fallback:land_binary_mask', 0)

        # default values for forcing variables
        o.set_config('environment:fallback:x_sea_water_velocity', 1.0) #m/s 
        o.set_config('environment:fallback:y_sea_water_velocity', 0) #m/s 
        o.set_config('environment:fallback:x_wind', 0 )#m/s
        o.set_config('environment:fallback:y_wind', 0) #m/s
        o.set_config('environment:fallback:sea_water_temperature', 12.0 )# degrees Celsius
        o.set_config('environment:fallback:sea_water_salinity', 35.0 )   # ppt
        o.set_config('environment:fallback:ocean_vertical_diffusivity', 0.000) # m2/s

        # use the update_terminal_velocity_constant() for the vertical motion, with daytime and nighttime positions
        o.set_config('biology:vertical_migration_speed_constant',1.0e-4) #'float(min=0.0, max=1e-3, default=None)', comment=' Constant vertical migration rate (m/s), if None, use values from update_terminal_velocity_pelagicegg()') #
        o.set_config('biology:vertical_position_daytime', -5.0)#'float(min=-1000.0, max=0.0, default=-5.0)',   comment='the depth a species is expected to inhabit during the day time, in meters, negative down') #
        o.set_config('biology:vertical_position_nighttime', -1.0) #'float(min=-1000.0, max=0.0, default=-1.0)', comment='the depth a species is expected to inhabit during the night time, in meters, negative down') #
        
        # set the "liveable" salinity and temperature thresholds so that constant values are within ranges.
        # the mortality will thus only be related to 'mortality_daily_rate'

        o.set_config('biology:temperature_min', 10.0)#'float(min=0.0, max=100.0, default=None)', comment=' lower threshold temperature where a species population quickly declines to extinction in degrees Celsius') #
        o.set_config('biology:temperature_max', 14.0)#'float(min=0.0, max=100.0, default=None)', comment=' upper threshold temperature where a species population quickly declines to extinction in degrees Celsius') #

        o.set_config('biology:salinity_min', 30.0)#'float(min=0.0, max=100.0, default=None)', comment=' lower threshold temperature where a species population quickly declines to extinction in degrees Celsius') #
        o.set_config('biology:salinity_max', 40.0)#'float(min=0.0, max=100.0, default=None)', comment=' upper threshold temperature where a species population quickly declines to extinction in degrees Celsius') #
    
        o.set_config('biology:mortality_daily_rate', 0.05)    # 'float(min=0.0, max=100.0, default=0.05)', comment='Mortality rate (percentage of biomass dying per day)') 

        o.seed_elements(lon=2, lat=60,z = -1.0, time=datetime(2021, 2, 24, 3, 46, 8), number=1)
        o.run(duration=timedelta(hours=4))
        
        self.assertAlmostEqual(o.elements.lon.max(), 2.258064243, 3)    
        self.assertAlmostEqual(o.elements.survival.max(), 0.9916926, 3)  # note: different vertical position at end of the simulation that previous test

    def test_PelagicPlankton_temperature_decay(self):

        o = PelagicPlanktonDrift(loglevel=0)
        #s.list_configspec()
        o.set_config('environment:fallback:land_binary_mask', 0)

        # default values for forcing variables
        o.set_config('environment:fallback:x_sea_water_velocity', 1.0) #m/s 
        o.set_config('environment:fallback:y_sea_water_velocity', 0) #m/s 
        o.set_config('environment:fallback:x_wind', 0 )#m/s
        o.set_config('environment:fallback:y_wind', 0) #m/s
        o.set_config('environment:fallback:sea_water_temperature', 13.0 )# degrees Celsius
        o.set_config('environment:fallback:sea_water_salinity', 35.0 )   # ppt
        o.set_config('environment:fallback:ocean_vertical_diffusivity', 0.000) # m2/s

        # use the update_terminal_velocity_constant() for the vertical motion, with daytime and nighttime positions
        o.set_config('biology:vertical_migration_speed_constant',1.0e-4) #'float(min=0.0, max=1e-3, default=None)', comment=' Constant vertical migration rate (m/s), if None, use values from update_terminal_velocity_pelagicegg()') #
        o.set_config('biology:vertical_position_daytime', -5.0)#'float(min=-1000.0, max=0.0, default=-5.0)',   comment='the depth a species is expected to inhabit during the day time, in meters, negative down') #
        o.set_config('biology:vertical_position_nighttime', -1.0) #'float(min=-1000.0, max=0.0, default=-1.0)', comment='the depth a species is expected to inhabit during the night time, in meters, negative down') #
        
        # set the "liveable" salinity and temperature thresholds so that constant values are within ranges.
        # the mortality will thus only be related to 'mortality_daily_rate'

        o.set_config('biology:temperature_min', 12.0)#'float(min=0.0, max=100.0, default=None)', comment=' lower threshold temperature where a species population quickly declines to extinction in degrees Celsius') #
        o.set_config('biology:temperature_max', 14.0)#'float(min=0.0, max=100.0, default=None)', comment=' upper threshold temperature where a species population quickly declines to extinction in degrees Celsius') #

        o.set_config('biology:salinity_min', 30.0)#'float(min=0.0, max=100.0, default=None)', comment=' lower threshold temperature where a species population quickly declines to extinction in degrees Celsius') #
        o.set_config('biology:salinity_max', 40.0)#'float(min=0.0, max=100.0, default=None)', comment=' upper threshold temperature where a species population quickly declines to extinction in degrees Celsius') #
    
        o.set_config('biology:mortality_daily_rate', 0.05)    # 'float(min=0.0, max=100.0, default=0.05)', comment='Mortality rate (percentage of biomass dying per day)') 

        o.seed_elements(lon=2, lat=60,z = -1.0, time=datetime(2021, 2, 24, 3, 46, 8), number=1)
        o.run(duration=timedelta(hours=4))
        # import pdb;pdb.set_trace()
        self.assertAlmostEqual(o.elements.lon.max(), 2.258064243, 3)    
        self.assertAlmostEqual(o.elements.survival.max(), 0.9916926, 3)  # note: different vertical position at end of the simulation that previous test
        length_run_in_days = 4./24.
        mortality_daily_rate = 0.05
        self.assertAlmostEqual(o.elements.survival.max(), 1-(length_run_in_days*mortality_daily_rate), 4)  # 1 - (length_run_in_days*mortality_daily_rate)

if __name__ == '__main__':
    unittest.main()
