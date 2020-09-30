#!/usr/bin/env python
"""
Generic plankton class
=============
"""

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.models.pelagicplankton_moana import PelagicPlanktonDrift
from opendrift.models.pelagicegg import PelagicEggDrift
from datetime import datetime
import numpy as np

o = PelagicPlanktonDrift(loglevel=20)  # Set loglevel to 0 for debug information
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
o.fallback_values['x_sea_water_velocity'] = 0 #m/s
o.fallback_values['y_sea_water_velocity'] = 0 #m/s 
o.fallback_values['x_wind'] = 0 #m/s
o.fallback_values['y_wind'] = 0 #m/s
o.fallback_values['sea_water_temperature'] = 12.0 # degrees Celsius
o.fallback_values['sea_water_salinity'] = 35.0    # ppt
o.fallback_values['ocean_vertical_diffusivity'] = 0.001 # m2/s
####################################################################################################################

####################################################################################################################
# SEEDING
# spawn NEA cod eggs at defined position and time
# time = datetime.utcnow()
o.seed_elements(lon=12.0, lat=68.3, z=np.linspace(0, -150, 100),
                radius=100, number=100,time=nordic_native.start_time, 
                diameter=0.0014, neutral_buoyancy_salinity=31.25) 

####################################################################################################################
# CONFIG
#%%
# Adjusting some configuration
o.set_config('general:coastline_action', 'previous')
o.set_config('drift:vertical_mixing', True)
#o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983') # windspeed parameterization for eddy diffusivity
o.set_config('vertical_mixing:diffusivitymodel', 'environment') # use eddy diffusivity from ocean model, or fallback value
o.set_config('vertical_mixing:timestep', 60.) # seconds - # Vertical mixing requires fast time step

if False:

    # plankton-specific config
    set_config('biology:constantIngestion', 0.5) #'float(min=0.0, max=1.0, default=0.5)', comment='Ingestion constant')
    set_config('biology:activemetabOn', 1.0) #'float(min=0.0, max=1.0, default=1.0)', comment='Active metabolism')
    set_config('biology:attenuationCoefficient', 0.18) #'float(min=0.0, max=1.0, default=0.18)', comment='Attenuation coefficient')
    set_config('biology:fractionOfTimestepSwimming',0.15 ) #'float(min=0.0, max=1.0, default=0.15)', comment='Fraction of timestep swimming')
    set_config('biology:lowerStomachLim', 0.3) #'float(min=0.0, max=1.0, default=0.3)', comment='Limit of stomach fullness for larvae to go down if light increases')
    set_config('biology:haddock', False) #'boolean(default=False)', comment='Species=haddock')
    set_config('biology:cod', True ) #'boolean(default=True)', comment='Species=cod')
####################################################################################################################

#%%
# Running model
o.run(steps=96, time_step=3600)

#%%
# Print and plot results.
# At the end the wind vanishes, and eggs come to surface
print(o)

o.plot(fast=True)
o.animation(fast=True, color='z')

#%%
# .. image:: /gallery/animations/example_codegg_0.gif

#%% Interactive slider (not working in browser)
o.plot_vertical_distribution()
