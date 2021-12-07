#!/usr/bin/env python

import os
import sys
import numpy as np
import shapefile
from datetime import datetime, timedelta
# from opendrift.readers import reader_ROMS_native_MOANA
from opendrift.models.fishlarvaeorient import FishLarvaeOrient

import time
start_time = time.time()

#direct, rheotaxis, cardinal, continuous_1, continuous_2
orient_option = sys.argv[1] # can be  'direct', 'rheotaxis', 'cardinal', 'continuous_1', 'continuous_2'

# 'direct'          : swims towards reef if close enough
# 'rheotaxis'       : swims against current
# 'cardinal'        : swims towards a specific cardinal heading 
# 'continuous_1'    : swims against current, then swims towards reef if close enough
# 'continuous_2'    : swims towards a specific cardinal heading   then swims towards reef if close enough

# model settings chosen so that in all cases, the larvae do feel the influence of user-defined habitat
# run times long enough so that particles go through different life stages

###############################
# MODEL SELECTION
###############################
o = FishLarvaeOrient(loglevel=0)
o.max_speed = 3
###############################
# READERS
###############################
o.set_config('general:use_auto_landmask', True) # prevent opendrift from making a new dynamical landmask
o.set_config('environment:fallback:x_sea_water_velocity', -.1) #  
o.set_config('environment:fallback:y_sea_water_velocity', -.1) #
###############################
# PARTICLE SEEDING
###############################
# 174.949513895 -40.029317969 28045 Taranaki Bight
o.seed_elements(174.949513895-.1, -40.029317969-.1, number=500, z=-2, time = [datetime(2021,1,1),datetime(2021,1,1)+timedelta(days=2.0)])
o.add_settlement_habitat('./habitat_moana/poly_habitat_taranaki.txt') # Location of the shapefile with the habitat
###############################
# PHYSICS
###############################
# horizontal and vertical diffusion
Kxy = 0.0 # m2/s-1
Kz = 0.001 # m2/s-1
o.set_config('environment:fallback:ocean_vertical_diffusivity', Kz) # specify constant ocean_vertical_diffusivity in m2.s-1
o.set_config('vertical_mixing:diffusivitymodel', 'constant') # constant >> use fall back values (can be environment (i.e. from reader), windspeed_Large1994 ,windspeed_Sundby1983)
o.set_config('drift:horizontal_diffusivity',Kxy) # using new config rather than current uncertainty
# seed
o.set_config('seed:ocean_only', True) # keep only particles from the "frame" that are on the ocean
# drift
o.set_config('drift:advection_scheme','runge-kutta4') # or 'runge-kutta'
o.set_config('drift:vertical_advection', True)
o.set_config('drift:vertical_mixing', True) 
o.set_config('vertical_mixing:timestep', 90.0)  # if some ocean_vertical_diffusivity!=0, turbulentmixing:timestep should be << 900 seconds
###############################
# Type of settlement
###############################
o.set_config('biology:settlement_in_habitat', True) # settlement restricted to suitable habitat only, specified by shapefile or txt file with nan-delimited polys
o.set_config('biology:min_settlement_age_seconds', 5*24*3600) # Beginning of the competency period when particle can settle
o.set_config('drift:max_age_seconds', 30*24*3600) # 
# coast/seafloor interactions set within module, to previous
###############################
# Vertical swimming
###############################
# vertical swimming speed and larval development stages
# enable the OVM for the larvae (need to turn on the vertical mixing (either constant or environment))
o.set_config('biology:OVM', True)                 # Ontogenetic Vertical Migration
o.set_config('biology:vertical_migration_speed_constant', 0.01) # vertical swimming speed in meters/s towards preferred depth
o.set_config('biology:pre_flexion',  0.5*24*3600)  # Beginning of the pre-flexion stage
o.set_config('biology:flexion',      2.0*24*3600)  # Beginning of the flexion stage. Beginning of the orientation 
o.set_config('biology:post_flexion', 4.0*24*3600)  # Beginning of the post-flexion stage
# Preferred depth per development stage
o.set_config('biology:depth_early_stage', -10.0)  # Preferred depth in meters before pre-flexion stage
o.set_config('biology:depth_pre_flexion', -15.0)  # Preferred depth in meters during pre-flexion stage
o.set_config('biology:depth_flexion', -20.0)      # Preferred depth in meters during flexion stage
o.set_config('biology:depth_post_flexion', -30.0) # Preferred depth in meters after flexion stage
o.set_config('biology:maximum_depth', -100.0)       # maximum depth of dispersal in meters (negative down)
###############################
# Orientation
###############################
# Larval orientation toward the nearest habitat
o.set_config('biology:orientation',orient_option) # orientation of the larvae: direct, rheotaxis, cardinal, continuous_1, continuous_2, or none
o.set_config('biology:beginning_orientation', 1.0*24*3600)   # Beginning of the flexion stage. Beginning of the orientation => if using OVM, values must agree
o.set_config('biology:max_orient_distance', 50.0)            # Orientation distance in kilometers. Maximum distance at which the reef habitat is detected - set very large to ensure influence
o.set_config('biology:cardinal_heading', 0.0)             # Cardinal heading for larval orientation when 'cardinal' option is selected
# Horizontal swimming speed:
o.set_config('biology:hatch_swimming_speed', 2.0)           # horizontal swimming speed at hatching, in cm/s
o.set_config('biology:settle_swimming_speed', 30.0)         # horizontal swimming speed at settlement, in cm/s

o.list_config()
# o.list_configspec()

###############################
# RUN 
###############################
# Running model (until end of driver data)
o.run(stop_on_error=False,
      time_step=timedelta(seconds = 900), # requires a small time-step to compute the orientation
      end_time = datetime(2021,1,1) + timedelta(days = 10.0),
      time_step_output=timedelta(seconds = 3600 * 3),
      outfile= 'example_fish_orient_%s.nc' % orient_option)

print(o)
if False:
      o.animation(fast=True, color='age_seconds') #corners=[163, 180, -52, -31])
      o.animation(fast=True, color='status') #corners=[163, 180, -52, -31])
      o.animation_profile()

o.animation(fast=True, color='age_seconds') #corners=[163, 180, -52, -31])
import matplotlib.pyplot as plt
plt.ion()
o.plot(fast=True)#, corners=[163, 180, -52, -31])
# add habitat polygons
import cartopy.crs as ccrs
gcrs = ccrs.PlateCarree()
plt.plot(o.multiShp[0].exterior.xy[0][:],o.multiShp[0].exterior.xy[1][:], transform = gcrs)
plt.plot(o.multiShp[1].exterior.xy[0][:],o.multiShp[1].exterior.xy[1][:], transform = gcrs)
import pdb;pdb.set_trace()
