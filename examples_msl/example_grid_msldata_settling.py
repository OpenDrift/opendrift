#!/usr/bin/env python

#  This script is based on https://github.com/OpenDrift/opendrift/tree/master/examples/example_grid.py
#  The aim is to test the use of reader_netCDF_MetOcean which allows using MetOcean data in OpenDrift
#   
#  This is also used as a reference case to compare Opendrift trajectories with ERcore trajectories
# 

import numpy as np
from datetime import datetime, timedelta
from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_netCDF_MetOcean
from opendrift.models.sedimentdrift3D import SedimentDrift3D
from opendrift.models.oceandrift3D import OceanDrift3D

###############################
# MODEL SELECTION
###############################
o = SedimentDrift3D(loglevel=0)  # Set loglevel to 0 for debug information
# o = OceanDrift3D(loglevel=0)  # Set loglevel to 0 for debug information

###############################
# READERS
###############################
reader_roms_cnz_depth = reader_netCDF_MetOcean.Reader('C:\metocean\cnz19800801_00z_surf.nc',variables_to_use = ['dep']) # i.e. only use the 'dep' variable in that file
reader_roms_cnz_depth.always_valid = True # this is useful if the time vector associated with the bathy does not fit with the simulation time 
#use always_valid for variables that can be used at all times, such as a bathy 

reader_roms_cnz = reader_netCDF_MetOcean.Reader('C:\metocean\cnz_surf_res200401.nc') # only uso,vso in this one - seems to be recognized ok 


# Making customised landmask (Basemap)

# reader_basemap = reader_basemap_landmask.Reader(
#                     llcrnrlon=3.5, llcrnrlat=59.9,
#                     urcrnrlon=5.5, urcrnrlat=61.2,
#                     resolution='h', projection='merc')

# o.add_reader([reader_roms_cnz,reader_roms_cnz_depth]) # ** need to add a reader with depth info when using settling
o.add_reader([reader_roms_cnz,reader_roms_cnz_depth]) # ** need to add a reader with depth info when using settling


 # no vertical diffusion infos available from readers : set fall_back constant values 
o.fallback_values['ocean_vertical_diffusivity'] = 0.0001

# all required variables that can be set using o.fall_back are generally listed 
# below the model class definition  e.g. see /models/sedimentdrift3D.py, line 36

###############################
# PARTICLE SEEDING
###############################


# Seeding some particles
# lons = np.linspace(3.5, 5.0, 100)
# lats = np.linspace(60, 61, 100)
# 
#  Point release
lon = 174.5133; lat = -41.2348; 
#  Rectangle release
# lons = np.linspace(170.0,170.5, 100) 
# lats = np.linspace(-39.5,-39.0, 100)
# lons, lats = np.meshgrid(lons, lats)
# lons = lons.ravel()
# lats = lats.ravel()
#

# Seed oil elements at defined position and time

# o.seed_elements(lons, lats, radius=0, number=10000,
#                 time=reader_roms_cnz.start_time)

o.seed_elements(lon, lat, radius=0, number=1000,time=datetime(2004,1,1),
                 z=0.0, terminal_velocity = -0.1) #, wind_drift_factor = 0, age_seconds = 0,)

# specific element variable such as terminal_velocity, can be specified here. 
# terminal_velocity>0 particle moves up, terminal_velocity<0 particle moves down
# 
# General variables to be defined lon,lat,z,time ..look in each /elements/element*.py for more properties
# 
# 

###############################
# PHYSICS
###############################

o.fallback_values['ocean_vertical_diffusivity'] = 0.0001 # specify constant ocean_vertical_diffusivity in m2.s-1

o.set_config('drift:current_uncertainty', 0.0)
o.set_config('drift:wind_uncertainty', 0.0)

o.set_config('processes:verticaladvection' , False) # no vertical current available, so no vertical advection
o.self.get_config('processes:resuspension',False) # already False be default but just for reference 
o.set_config('processes:turbulentmixing', True) # 
o.set_config('turbulentmixing:diffusivitymodel', 'environment') # i.e. specified from model or constant
o.set_config('turbulentmixing:TSprofiles',False)
o.set_config('turbulentmixing:timestep', 1800) 
# A potential workaround to use a given settling velocity : terminal_velocity, and no added vertical diffusion is to use :
# o.fallback_values['ocean_vertical_diffusivity'] = 0.0
# The time step governing the vertical settling (and mixing if ocean_vertical_diffusivity~=0) is :
# o.set_config('turbulentmixing:timestep', 1800)  # can be set to same as "run" time step to reproduce ERcore behaviour


###############################
# RUN 
###############################

# Running model (until end of driver data)
o.run(time_step=1800, end_time = datetime(2004,1,2), outfile='opendrift_settling_withdepth.nc',time_step_output = 1800)

# the start time is defined by seed_elements, the end_time is defined by either steps=number of step, duration = timedelta, or end_time= datetime

###############################
# PLOTS / ANIMATION
###############################

# Print and plot results
print o

o.animation()
o.plot()
