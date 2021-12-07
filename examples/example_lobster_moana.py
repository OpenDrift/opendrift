#!/usr/bin/env python

import os
import sys
import numpy as np
import shapefile
from datetime import datetime, timedelta
# from opendrift.readers import reader_ROMS_native_MOANA
from opendrift.models.lobsterlarvae import LobsterLarvae

###############################
# MODEL SELECTION
###############################
o = LobsterLarvae(loglevel=0)
o.max_speed = 3
###############################
# READERS
###############################
o.set_config('general:use_auto_landmask', True)
# constant current
o.set_config('environment:fallback:x_sea_water_velocity', -.1) #  
o.set_config('environment:fallback:y_sea_water_velocity', -.1) #
###############################
# PARTICLE SEEDING
###############################
o.seed_elements(174.949513895-.1, -40.029317969-.1, number=500, z=-2, time = [datetime(2021,1,1),datetime(2021,1,1)+timedelta(days=2.0)]) #
o.add_settlement_habitat('./habitat_moana/poly_habitat_taranaki.txt') # Location of the shapefile or text file with the habitat polygon(s)
###############################
# PHYSICS
###############################
# horizontal and vertical diffusion
Kxy = 0.0 # m2/s-1
Kz = 0.001 # m2/s-1
o.set_config('environment:fallback:ocean_vertical_diffusivity', Kz) # specify constant ocean_vertical_diffusivity in m2.s-1
o.set_config('drift:advection_scheme','runge-kutta4') # or 'runge-kutta'
o.set_config('drift:vertical_mixing', True) 
o.set_config('vertical_mixing:diffusivitymodel', 'constant') # constant >> use fall back values
o.set_config('vertical_mixing:timestep', 90.0)  
o.set_config('vertical_mixing:TSprofiles',False)
o.set_config('drift:horizontal_diffusivity',Kxy) 
###############################
# Life stage details
###############################
o.set_config('drift:max_age_seconds', 30*24*3600) # maximum possible age of particles
# Orientation
# Larval orientation toward the nearest habitat
o.set_config('biology:settlement_in_habitat', True)             # settlement restricted to suitable habitat only, specified by shapefile or txt file with nan-delimited polys
o.set_config('biology:direct_orientation_habitat',True)         # particle will swim towards its habitat when age>age_beginning_orientation
o.set_config('biology:age_beginning_orientation', 1.0*24*3600)  # Beginning of the orientation stage...could be linked to phyllosoma stages if relevant 
o.set_config('biology:max_orient_distance', 50.0)               # Maximum distance at which the reef habitat is detected - set very large to ensure influence in example
o.set_config('biology:min_settlement_age_seconds', 1*24*3600.0) # Minimum age before beaching/settling can occur in user-define habitat (from bivalve module)
# life stage
o.set_config('biology:min_swimming_speed_puerulus',10.0) # in cm/s minimum swimming speed of the puerulus when cruising toward the habitat
o.set_config('biology:max_swimming_speed_puerulus',30.0) # in cm/s maximum swimming speed of the puerulus when cruising toward the habitat
o.set_config('biology:mid_stage_phyllosoma', 2*24*3600)  # start age for 'middle stage phyllosoma'
o.set_config('biology:late_stage_phyllosoma',4*24*3600)  # start age for 'late stage phyllosoma'
# vertical motion
o.set_config('biology:vertical_position_daytime', -5.0)        # the depth a species is expected to inhabit during the day time, in meters, negative down') #
o.set_config('biology:vertical_position_nighttime', -1.0)      # the depth a species is expected to inhabit during the night time, in meters, negative down') #
o.set_config('biology:vertical_migration_speed_constant',1e-3) #'Constant vertical migration rate (m/s)
o.set_config('biology:maximum_larvae_depth', -100.0)           # maximum depth of dispersal in meters (negative down)

o.list_config()
# o.list_configspec()

###############################
# RUN 
###############################
# Running model (until end of driver data)
o.run(stop_on_error=False,
      time_step=timedelta(seconds = 900), # requires a small time-step to compute the orientation
      end_time = datetime(2021,1,1) + timedelta(days = 10.0),
      time_step_output=timedelta(seconds = 3600 * 3))
      # outfile= 'example_lobster_orient.nc')

print(o)

o.animation(fast=True, color='z')
import matplotlib.pyplot as plt
plt.ion()
o.plot(fast=True)#, corners=[163, 180, -52, -31])
# add habitat polygons
import cartopy.crs as ccrs
gcrs = ccrs.PlateCarree()
plt.plot(o.multiShp[0].exterior.xy[0][:],o.multiShp[0].exterior.xy[1][:], transform = gcrs)
plt.plot(o.multiShp[1].exterior.xy[0][:],o.multiShp[1].exterior.xy[1][:], transform = gcrs)
import pdb;pdb.set_trace()
