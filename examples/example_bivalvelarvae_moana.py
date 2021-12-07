#!/usr/bin/env python
"""
Bivalve Larvae class
=============
"""


from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_ROMS_native_MOANA
from datetime import datetime,timedelta
from opendrift.models.bivalvelarvae import BivalveLarvae
import numpy as np
import sys

# call the script using 
# python example_bivalvelarvae_moana.py 0, for example with no habitat
# python example_bivalvelarvae_moana.py 1  for example with habitat

use_habitat = int(sys.argv[1])

o = BivalveLarvae(loglevel=0)  # Set loglevel to 0 for debug information

###############################
# READERS
###############################
thredds_path_1 = 'http://thredds.moanaproject.org:8080/thredds/dodsC/moana/ocean/NZB/v1.9/monthly_avg/nz5km_avg_201603.nc'
reader_moana_v19 = reader_ROMS_native_MOANA.Reader(thredds_path_1) # load data for that year
reader_moana_v19.multiprocessing_fail = True # bypass the use of multi core for coordinates conversion and seems to make the model run much faster.
# use native landmask of ROMS files
o.add_reader([reader_moana_v19]) # [reader_landmask,reader_moana_v19]
o.set_config('general:use_auto_landmask', False) # prevent opendrift from making a new dynamical landmask
###############################
# PARTICLE SEEDING
###############################
# 174.949513895 -40.029317969 28045 Taranaki Bight
o.seed_elements( lon = 174.949513895-.1, lat = -40.029317969-.1, 
                 number=100, z = - np.linspace(2,50,100) ,
                 radius = 1000,
                 terminal_velocity = -1.0, # setting a large settling velocity on purpose to make seafloor contacts
                 time = reader_moana_v19.start_time)
###############################
# PHYSICS & CONFIG
###############################
# horizontal and vertical diffusion
Kxy = 0.1 # m2/s-1
Kz = 0.001 # m2/s-1 - setting a large 
o.set_config('environment:fallback:ocean_vertical_diffusivity', Kz) # specify constant ocean_vertical_diffusivity in m2.s-1
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'constant') # constant >> use fall back values (can be environment (i.e. from reader), windspeed_Large1994 ,windspeed_Sundby1983)
o.set_config('vertical_mixing:timestep', 900.) # seconds - # Vertical mixing requires fast time step  (but for constant diffusivity, use same as model step)
o.set_config('drift:horizontal_diffusivity',Kxy) # using new config rather than current uncertainty
o.set_config('drift:max_age_seconds', 10*24*3600) # max age before dying
###############################
# bivalve specific
###############################
o.set_config('biology:min_settlement_age_seconds', 2*24*3600) # minimum age before settling can occur, on seafloor, coast or habitat

if use_habitat ==1 :
    o.set_config('biology:settlement_in_habitat', True) # settlement restricted to suitable habitat only, specified by shapefile or txt file with nan-delimited polys
    o.add_settlement_habitat('./habitat_moana/poly_habitat_taranaki_bivalve.txt') # Location of the files with the user-defined habitat polygons
                                                                                  # habitat files can be text file with simple nan-delimited polygon(s) 
else:
    o.set_config('biology:settlement_in_habitat', False) # particle will settle on coast or seafloor after it reached 'min_settlement_age_seconds'

o.list_config()

###############################
# RUN 
###############################
# Running model (until end of driver data)
o.run(stop_on_error=False,
      time_step=timedelta(seconds = 900),
      end_time = reader_moana_v19.start_time + timedelta(days = 4.0),
      time_step_output=timedelta(seconds = 3600 * 3))

import matplotlib.pyplot as plt
plt.ion()
o.plot(fast=True)#, corners=[163, 180, -52, -31])
if use_habitat ==1 :
    # add habitat polygons
    import cartopy.crs as ccrs
    gcrs = ccrs.PlateCarree()
    plt.plot(o.multiShp.exterior.xy[0][:],o.multiShp.exterior.xy[1][:], transform = gcrs)

import pdb;pdb.set_trace()
# o.animation(fast=True, color='status') #corners=[163, 180, -52, -31])
# o.animation_profile()
# o.plot(fast=True)
