#!/usr/bin/env python

import os
import sys
import numpy as np
from datetime import datetime, timedelta
# from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_ROMS_native_MOANA
from opendrift.models.oceandrift import OceanDrift
# from opendrift.models.bivalvelarvae import BivalveLarvae

###############################
# MODEL SELECTION
###############################
o = OceanDrift(loglevel=0)
o.max_speed = 3.0#
###############################
# READERS
###############################

thredds_path = 'http://thredds.moanaproject.org:8080/thredds/dodsC/moana/ocean/NZB/v1.9/raw_3D/nz5km_his_201712.nc'
reader_moana_v19 = reader_ROMS_native_MOANA.Reader(thredds_path) # load data for that year
reader_moana_v19.multiprocessing_fail = False # bypass the use of multi core for coordinates conversion and seems to make the model run much faster.

# # Making customised landmask - not required here, using ROMS landmask 
# reader_landmask = reader_global_landmask.Reader(
#                     llcrnrlon=171.0, llcrnrlat=184.5,
#                     urcrnrlon=-42.0, urcrnrlat=-34.0)
# max is 185deg

# use native landmask of ROMS files
o.add_reader([reader_moana_v19]) # [reader_landmask,reader_moana_v19]
o.set_config('general:use_auto_landmask', False) # prevent opendrift from making a new dynamical landmask

###############################
# PARTICLE SEEDING
###############################
# generic point release location for test
lon0 = 177.269014281   
lat0 = -37.8719 
nb_parts = 1000
radius_in_m = 1000

# frame release
# lons,lats = np.meshgrid(np.linspace(173.5, 173.8, 10),np.linspace(-42.6,-42.3, 10))
# lon0 = np.ravel(lons)
# lat0 = np.ravel(lats)
# nb_parts = 100
# radius_in_m = 0

z = np.random.uniform(-50,0,size=nb_parts) # generate random depth
# continuous release from tstart_release to tend_release within a 1000m radius
o.seed_elements(lon0, lat0, 
                      radius=radius_in_m, 
                      number=nb_parts, 
                      z=z,
                      time=[reader_moana_v19.start_time, reader_moana_v19.start_time + timedelta(days=1.0)]) #

###############################
# PHYSICS
###############################
o.set_config('environment:fallback:x_wind', 0.0)
o.set_config('environment:fallback:y_wind', 0.0)
o.set_config('environment:fallback:x_sea_water_velocity', 0.0)
o.set_config('environment:fallback:y_sea_water_velocity', 0.0)
o.set_config('environment:fallback:sea_floor_depth_below_sea_level', 100000.0)

# horizontal and vertical diffusion
Kxy = 1.0 #0.1176 # m2/s-1
Kz = 0.001 # m2/s-1
o.set_config('environment:fallback:ocean_vertical_diffusivity', Kz) # specify constant ocean_vertical_diffusivity in m2.s-1
o.set_config('vertical_mixing:diffusivitymodel', 'constant') # constant >> use fall back values o.set_config('environment:fallback:ocean_vertical_diffusivity'] = Kz for all profile
# can be environment (i.e. from reader), or  windspeed_Large1994 ,windspeed_Sundby1983, 
o.set_config('drift:horizontal_diffusivity',Kxy) # using new config rather than current uncertainty
# seed
o.set_config('seed:ocean_only',True) # keep only particles from the "frame" that are on the ocean
# drift
o.set_config('drift:advection_scheme','runge-kutta4') # or 'runge-kutta'
o.set_config('drift:current_uncertainty', 0.0 ) # note current_uncertainty can be used to replicate an horizontal diffusion spd_uncertain = sqrt(Kxy*2/dt)  
o.set_config('drift:max_age_seconds', 10*24*3600) # 
o.set_config('drift:vertical_advection', False) 
o.set_config('drift:vertical_mixing', True) 
o.set_config('vertical_mixing:timestep', 90.0)  # if some ocean_vertical_diffusivity!=0, turbulentmixing:timestep should be << 900 seconds

o.set_config('general:coastline_action','stranding') # option('none', 'stranding', 'previous', default='stranding')
o.set_config('general:seafloor_action','lift_to_seafloor')

o.disable_vertical_motion()  #Deactivate any vertical processes/advection"""

o.list_config()
# o.list_configspec()

###############################
# RUN 
###############################
# Running model (until end of driver data)
o.run(stop_on_error=False,
      time_step=900, 
      end_time = reader_moana_v19.start_time+timedelta(days=3.0), 
      outfile= 'opendrift_using_moana.nc',
      time_step_output = 3600.0)

o.plot()