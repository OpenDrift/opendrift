#!/usr/bin/env python

#  This script is based on https://github.com/OpenDrift/opendrift/tree/master/examples/example_grid.py
#  The aim is to test the use of reader_netCDF_MetOcean which allows using MetOcean data in OpenDrift
#   
#  This is also used as a reference case to compare Opendrift trajectories with ERcore trajectories
# 
#
#  This example show how to use the SedimentDrift3D model, which allows modelling horizontal and vertical advection due to currents
#  horizontal and vertical diffusion (mixing), and buoyancy-related setlling 
#
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
data_pth = 'C:/Users/simon/Google Drive/R&D/OpenDrift/data_testing/'

reader_ww3_wave = reader_netCDF_MetOcean.Reader(data_pth + 'ww3_indiano_20150101_201501007.nc') # Use all variables 
reader_cfsr2_current = reader_netCDF_MetOcean.Reader(data_pth + 'cfsr2_uv_20150101_201501007.nc')  # residual surface current only

# Making customised landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=20.0, llcrnrlat=-60.0,
                    urcrnrlon=60.0, urcrnrlat=-0.0,
                    resolution='l', projection='merc') # resolution can be c (crude, the default), l (low), i (intermediate), h (high), f (full)

o.add_reader([reader_basemap,reader_cfsr2_current]) # reader_ww3_wave

###############################
# PARTICLE SEEDING
###############################

# Seeding some particles
# lons = np.linspace(3.5, 5.0, 100)
# lats = np.linspace(60, 61, 100)
# 
#  Point release
lon = 39.1915
lat = -44.3692

o.seed_elements(lon, lat, radius=0, number=10,time=reader_cfsr2_current.start_time,
                 z=0.0, terminal_velocity = -0.000) # wind_drift_factor = 0.04 ) # age_seconds = 0,)

# specific element variable such as terminal_velocity, can be specified here. 
# terminal_velocity>0 particle moves up, terminal_velocity<0 particle moves down
# 
# General variables to be defined lon,lat,z,time ..look in each /elements/element*.py for more properties
# 
# 

###############################
# PHYSICS
###############################

# these will list all possible options for that model


# diffusion - constant in that example
o.fallback_values['ocean_horizontal_diffusivity'] = 0.0 # specify constant ocean_horizontal_diffusivity in m2.s-1
o.fallback_values['ocean_vertical_diffusivity'] = 0.0000 # specify constant ocean_vertical_diffusivity in m2.s-1

# drift
o.set_config('drift:scheme','euler') # or 'runge-kutta'
o.set_config('drift:current_uncertainty', 0.0)
o.set_config('drift:wind_uncertainty', 0.0)
o.set_config('drift:stokes_drift', True)
o.set_config('drift:use_tabularised_stokes_drift', True)
o.set_config('drift:tabularised_stokes_drift_fetch','25000')

#processes
o.set_config('processes:verticaladvection' , False) # no vertical current available, so no vertical advection
o.set_config('processes:resuspension',False) # already False be default but just for reference 
o.set_config('processes:turbulentmixing', False) 
# Note that turbulentmixing include the buoyancy-related settling for now (terminal_velocity)
# if False, no vertical mixing NOR settling will occur - might need to be changed to allow pure buoyancy-drive settling and skip the vertical mixing computations?
# 
# 
o.set_config('turbulentmixing:diffusivitymodel', 'environment') # i.e. specified from model or constant
o.set_config('turbulentmixing:TSprofiles',False)
o.set_config('turbulentmixing:timestep', 1800)  


o.list_config()
o.list_configspec()

###############################
# RUN 
###############################
# import pdb;pdb.set_trace()
# Running model (until end of driver data)
o.run(time_step=1800, end_time = reader_cfsr2_current.start_time + timedelta(days = 7.0), outfile='opendrift_sedimentdrift3d_all_inputs.nc',time_step_output = 1800)
# the start time is defined by seed_elements, the end_time is defined by either steps=number of step, duration = timedelta, or end_time= datetime

###############################
# PLOTS / ANIMATION
###############################


# Print and plot results
print o
# o.animation()
o.plot()
