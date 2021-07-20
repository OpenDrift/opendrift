#!/usr/bin/env python

import os
import sys
import numpy as np
import shapefile
from datetime import datetime, timedelta
from opendrift.readers import reader_ROMS_native_MOANA
#from opendrift.models.oceandrift import OceanDrift
from opendrift.models.bivalvelarvae import BivalveLarvae
import time
start_time = time.time()

####
# Month selection
####

runtime = [datetime(2016,3,2), datetime(2016,3,3)]
finish_time = (datetime(2016,3,2) + timedelta(days = 31.0))


###############################
# MODEL SELECTION
###############################
o = BivalveLarvae(loglevel=100)
o.max_speed = 3#
###############################
# READERS
###############################
thredds_path_1 = 'http://thredds.moanaproject.org:8080/thredds/dodsC/moana/ocean/NZB/v1.9/monthly_avg/nz5km_avg_201603.nc'
reader_moana_v19 = reader_ROMS_native_MOANA.Reader(thredds_path_1) # load data for that year
reader_moana_v19.multiprocessing_fail = True # bypass the use of multi core for coordinates conversion and seems to make the model run much faster.

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

# Define the starting position of the particles within the buffer of intertidal rocky shore on the West Coast of the North Island (NZ)
nb_parts = 100
points = np.loadtxt('./Release_centroid_nat_dist_paua.xyz', delimiter='\t', dtype=str)
plon = points[:,0].astype(np.float)
plat = points[:,1].astype(np.float)
tot_parts = nb_parts * len(plon) # total number of particles released
plon = np.tile(plon, nb_parts)
plat = np.tile(plat, nb_parts)

# Define the release time. author: Calvin Quigley (19/04/2021)
def create_seed_times(start, end, delta):
  """
  crate times at given interval to seed particles
  """
  out = []
  start_t = start
  end_t = datetime.strptime(str(end), "%Y-%m-%d %H:%M:%S")
  while start_t < end:
    out.append(start_t) 
    start_t += delta
  return out

times = create_seed_times(runtime[0], 
                          runtime[1], timedelta(days = 1))

# Define release depth
z = np.random.uniform(-0.5,-8,size=tot_parts) # generate random depth
#z = 'seafloor+1'

# Seed particles
for i in range(len(times)):
	o.seed_elements(plon, plat, number=tot_parts, z=z, time = times[i], terminal_velocity=0.001)

# seed
o.set_config('seed:ocean_only',False) # keep only particles from the "frame" that are on the ocean


#lon0 = 173.11611246985
#lat0 = -35.3876529758404
#radius_in_m = 1000
#nb_parts = 1000
#z = np.random.uniform(-50,0,size=nb_parts) # generate random depth
# continuous release from tstart_release to tend_release within a 1000m radius
#o.seed_elements(plon, plat, 
#                      #radius=radius_in_m, 
#                      number=tot_parts, 
#                      z=z,
#                      time=[reader_moana_v19.start_time, reader_moana_v19.start_time + timedelta(days=1.0)]) #

###############################
# PHYSICS
###############################
o.set_config('environment:fallback:x_wind', 0.0)
o.set_config('environment:fallback:y_wind', 0.0)
o.set_config('environment:fallback:x_sea_water_velocity', 0.0)
o.set_config('environment:fallback:y_sea_water_velocity', 0.0)
o.set_config('environment:fallback:sea_floor_depth_below_sea_level', 12000)

# horizontal and vertical diffusion
Kxy = 0.1176 #0.1176 # m2/s-1
Kz = 0.01 # m2/s-1
o.set_config('environment:fallback:ocean_vertical_diffusivity', Kz) # specify constant ocean_vertical_diffusivity in m2.s-1
o.set_config('vertical_mixing:diffusivitymodel', 'constant') # constant >> use fall back values o.set_config('environment:fallback:ocean_vertical_diffusivity'] = Kz for all profile
# can be environment (i.e. from reader), or  windspeed_Large1994 ,windspeed_Sundby1983, 
o.set_config('drift:horizontal_diffusivity',Kxy) # using new config rather than current uncertainty
# seed
o.set_config('seed:ocean_only', True) # keep only particles from the "frame" that are on the ocean
# drift
o.set_config('drift:advection_scheme','runge-kutta4') # or 'runge-kutta'
o.set_config('drift:current_uncertainty', 0.0 ) # note current_uncertainty can be used to replicate an horizontal diffusion spd_uncertain = sqrt(Kxy*2/dt)  
o.set_config('drift:vertical_advection', True)
o.set_config('vertical_mixing:timestep', 90.0)  # if some ocean_vertical_diffusivity!=0, turbulentmixing:timestep should be << 900 seconds
o.set_config('drift:vertical_mixing', True) 

###############################
# Type of settlement
###############################
o.habitat('./habitat/National_distribution_paua.shp') # Location of the shapefile with the habitat
o.set_config('drift:settlement_in_habitat', True)
o.set_config('drift:max_age_seconds', 9*24*3600) # 
o.set_config('drift:min_settlement_age_seconds', 5*24*3600)
o.set_config('general:seafloor_action', 'lift_to_seafloor')
o.set_config('drift:Haliotis_iris', True)

###############################
# Vertical swimming
###############################
# Need a null terminal_velocity
o.set_config('drift:active_vertical_swimming', False) # Correlated random walk across the water column when advected away from coastal habitats
#o.set_config('drift:vertical_velocity', (0.001/10)) # Vertical swimming speed of the larvae: in meter/seconds
o.set_config('drift:maximum_depth', -50) # Maximum depth of larvae: negative in meters
#o.set_config('drift:persistence', 50) # Control the persistence (memory) of the vertical movement to create a correlated random walk and sample the water column


o.list_config()
# o.list_configspec()

###############################
# RUN 
###############################
# Running model (until end of driver data)
o.run(stop_on_error=False,
      time_step=timedelta(seconds = 900), 
      end_time = finish_time,
      time_step_output=timedelta(minutes = 60*24),
      outfile= 'test_development_Opendrift_using_moana.nc')

print("--- %s seconds ---" % (time.time() - start_time))
print(process.memory_info().rss)  # in bytes 


print(o)
o.plot(fast=True, color='z', corners=[163, 180, -52, -31])
o.animation(fast=True, color='z', corners=[163, 180, -52, -31])
o.animation_profile()
