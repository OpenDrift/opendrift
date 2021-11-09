#!/usr/bin/env python
"""
Bivalve Larvae class
=============
"""

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from datetime import datetime,timedelta
from opendrift.models.bivalvelarvae import BivalveLarvae
import numpy as np

o = BivalveLarvae(loglevel=0)  # Set loglevel to 0 for debug information
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
o.set_config('environment:fallback:ocean_vertical_diffusivity',0.001) # m2/s
####################################################################################################################

####################################################################################################################
# SEEDING
# time = datetime.utcnow()
o.seed_elements(lon=12.0, lat=68.3, z=-3.0,
                radius=100, number=100,
                time=nordic_native.start_time,neutral_buoyancy_salinity=31.25,
                terminal_velocity = 0.0)  

####################################################################################################################
# CONFIG
# Adjusting some configuration
o.set_config('general:coastline_action', 'previous')
o.set_config('drift:vertical_advection', True)
o.set_config('drift:vertical_mixing', True)
o.set_config('biology:min_settlement_age_seconds', 3600) # minimum age before settling can occur

o.set_config('vertical_mixing:diffusivitymodel', 'constant') # use fallback value ocean_vertical_diffusivity
o.set_config('vertical_mixing:timestep', 900.) # seconds - # Vertical mixing requires fast time step  (but for constant diffusivity, use same as model step)
####################################################################################################################

# Running model
o.run(end_time=nordic_native.start_time + timedelta(days=4.0), time_step=900)

# Print and plot results.
# At the end the wind vanishes, and eggs come to surface
print(o)


o.plot(fast=True)

o.animation(fast=True, color='z')

# Interactive slider (not working in browser)
o.plot_vertical_distribution()
