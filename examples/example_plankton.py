#!/usr/bin/env python
"""
Generic plankton class
=============
"""

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.models.pelagicplankton_moana import PelagicPlanktonDrift
from opendrift.models.pelagicegg import PelagicEggDrift
from datetime import datetime,timedelta
import numpy as np

o = PelagicPlanktonDrift(loglevel=0)  # Set loglevel to 0 for debug information
# o = PelagicEggDrift(loglevel=20)  # Set loglevel to 0 for debug information

####################################################################################################################
# READERS
# 
# Atmospheric model for wind
# reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
#     '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Ocean model for current
# reader_norkyst = reader_netCDF_pytCF_generic.Reader(o.test_data_folder() +
#     '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
# # Creating and adding reader for Nordic 4km current dataset
nordic_native = reader_ROMS_native.Reader(o.test_data_folder() +
    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
o.add_reader(nordic_native)

# default values for forcing variables
o.set_config('environment:fallback:x_sea_water_velocity', 0) #m/s 
o.set_config('environment:fallback:y_sea_water_velocity', 0) #m/s 
o.set_config('environment:fallback:x_wind', 0 )#m/s
o.set_config('environment:fallback:y_wind', 0) #m/s
o.set_config('environment:fallback:sea_water_temperature', 12.0 )# degrees Celsius
o.set_config('environment:fallback:sea_water_salinity', 35.0 )   # ppt
o.set_config('environment:fallback:ocean_vertical_diffusivity', 0.001) # m2/s
####################################################################################################################

####################################################################################################################
# SEEDING
# spawn NEA cod eggs at defined position and time
# time = datetime.utcnow()
o.seed_elements(lon=12.0, lat=68.3, z=-3.0, #z should be between vertical_position_daytime and vertical_position_daytime z=np.linspace(0, -150, 100),
                radius=100, number=100,time=nordic_native.start_time, 
                diameter=0.0014, neutral_buoyancy_salinity=31.25)  

# o.seed_elements(lon=12.0, lat=68.3, z=np.linspace(0, -150, 100),
#                 radius=100, number=100,time=nordic_native.start_time, 
#                 diameter=0.0014, neutral_buoyancy_salinity=31.25) 

####################################################################################################################
# CONFIG
# Adjusting some configuration
o.set_config('general:coastline_action', 'previous')
o.set_config('drift:vertical_advection', True)

o.set_config('drift:vertical_mixing', True)
# o.o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983') # windspeed parameterization for eddy diffusivity
o.set_config('vertical_mixing:diffusivitymodel', 'environment') # use eddy diffusivity from ocean model, or fallback value
o.set_config('vertical_mixing:timestep', 3600.) # seconds - # Vertical mixing requires fast time step  (but for constant diffusivity, use same as model step)
####################################################################################################################
if True:
    # plankton-specific config for developped module based on 
    o.set_config('biology:mortality_daily_rate', 0.05)    # 'float(min=0.0, max=100.0, default=0.05)', comment='Mortality rate (percentage of biomass dying per day)') 
    o.set_config('biology:min_settlement_age_seconds', 5*24*3600.0)  #'float(min=0.0, max=100.0, default=0.0)', comment='Minimum age before beaching can occur, in seconds')
    o.set_config('biology:vertical_position_daytime', -5.0)#'float(min=-1000.0, max=0.0, default=-5.0)',   comment='the depth a species is expected to inhabit during the day time, in meters, negative down') #
    o.set_config('biology:vertical_position_nighttime', -1.0) #'float(min=-1000.0, max=0.0, default=-1.0)', comment='the depth a species is expected to inhabit during the night time, in meters, negative down') #
    o.set_config('biology:vertical_migration_speed_constant',1e-4) #'float(min=0.0, max=1e-3, default=None)', comment=' Constant vertical migration rate (m/s), if None, use values from update_terminal_velocity_pelagicegg()') #
    o.set_config('biology:temperature_min', 5.0)#'float(min=0.0, max=100.0, default=None)', comment=' lower threshold temperature where a species population quickly declines to extinction in degrees Celsius') #
    o.set_config('biology:temperature_max', 25.0)#'float(min=0.0, max=100.0, default=None)', comment=' upper threshold temperature where a species population quickly declines to extinction in degrees Celsius') #
    o.set_config('biology:temperature_tolerance', 1.0)#'float(min=0.0, max=1.0, default=1.0)', comment=' temperature tolerance before dying in degrees Celsius') #
    o.set_config('biology:salinity_min', 30.0)#'float(min=0.0, max=100.0, default=None)', comment=' lower threshold salinity where a species population quickly declines to extinction in ppt') #
    o.set_config('biology:salinity_max', 39.0)#'float(min=0.0, max=100.0, default=None)', comment=' upper threshold salinity where a species population quickly declines to extinction in ppt') #
    o.set_config('biology:salinity_tolerance',1.0)#'float(min=0.0, max=1.0, default=1.0)', comment=' salinity tolerance before dying in ppt') #

    # to switch off the constant migration rate towards day or night time position, and use update_terminal_velocity_pelagicegg() :
    # o.set_config('biology:vertical_migration_speed_constant',None) 
####################################################################################################################
if False:
    # in dev...
    # plankton-specific config for module in :
    # https://github.com/trondkr/KINO-ROMS/blob/master/Romagnoni-2019-OpenDrift/kino/pelagicplankton.py
    o.set_config('biology:constantIngestion', 0.5) #'float(min=0.0, max=1.0, default=0.5)', comment='Ingestion constant')
    o.set_config('biology:activemetabOn', 1.0) #'float(min=0.0, max=1.0, default=1.0)', comment='Active metabolism')
    o.set_config('biology:attenuationCoefficient', 0.18) #'float(min=0.0, max=1.0, default=0.18)', comment='Attenuation coefficient')
    o.set_config('biology:fractionOfTimestepSwimming',0.15 ) #'float(min=0.0, max=1.0, default=0.15)', comment='Fraction of timestep swimming')
    o.set_config('biology:lowerStomachLim', 0.3) #'float(min=0.0, max=1.0, default=0.3)', comment='Limit of stomach fullness for larvae to go down if light increases')
    o.set_config('biology:haddock', False) #'boolean(default=False)', comment='Species=haddock')
    o.set_config('biology:cod', True ) #'boolean(default=True)', comment='Species=cod')
####################################################################################################################


# Running model
o.run(end_time=nordic_native.start_time + timedelta(hours=4.0), time_step=3600)

# Print and plot results.
# At the end the wind vanishes, and eggs come to surface
print(o)

o.plot(fast=True)

o.animation(fast=True, color='z')

# Interactive slider (not working in browser)
o.plot_vertical_distribution()
