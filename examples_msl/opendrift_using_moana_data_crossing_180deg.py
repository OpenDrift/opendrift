#!/usr/bin/env python
# 
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
###############################
# READERS
###############################
thredds_path = 'http://thredds.moanaproject.org:8080/thredds/dodsC/moana/ocean/NZB/v1.9/raw_3D/nz5km_his_201712.nc'
reader_moana_v19 = reader_ROMS_native_MOANA.Reader(thredds_path) # load data for that year
reader_moana_v19.multiprocessing_fail = True # bypass the use of multi core for coordinates conversion and seems to make the model run much faster.
# use native landmask of ROMS files
o.add_reader([reader_moana_v19]) # [reader_landmask,reader_moana_v19]
o.set_config('general:use_auto_landmask', False) # prevent opendrift from making a new dynamical landmask, use ROMS landmask instead

###############################
# PARTICLE SEEDING
###############################
#  Rectangle release
lon = np.linspace(179.0,181.0, 100)
lat = np.linspace(-45.,-46., 10)
# these lon/lat spread across both roms_nz and cfsr data
depth = -np.linspace(0.,0., 1000) # surface=0, negative down
lons, lats = np.meshgrid(lon, lat)
lons = lons.ravel()
lats = lats.ravel()

o.seed_elements(lons.ravel(), lats.ravel(), 
                z = depth,
                time=reader_moana_v19.start_time) #

###############################
# PHYSICS
###############################
o.disable_vertical_motion()  #Deactivate any vertical processes/advection"""op
###############################
# RUN 
###############################
# Running model (until end of driver data)
o.run(stop_on_error=False,
      time_step=900, 
      end_time = reader_moana_v19.start_time+timedelta(days=3.0), 
      outfile= 'opendrift_using_moana_crossing180deg.nc',
      time_step_output = 3600.0)

# o.plot()