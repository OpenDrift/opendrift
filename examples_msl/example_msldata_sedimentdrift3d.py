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

reader_roms_cnz_depth = reader_netCDF_MetOcean.Reader(data_pth + 'current_roms_cnz_dep.nc',variables_to_use = ['dep']) # i.e. only use the 'dep' variable in that file
reader_roms_cnz_depth.always_valid = True # this is useful if the time vector associated with the bathy does not fit with the simulation time 
#use always_valid for variables that can be used at all times, such as a bathy 

reader_roms_cnz_dav = reader_netCDF_MetOcean.Reader(data_pth + 'current_roms_cnz_surf_20050101_20050107.nc',variables_to_use = ['dep','um','vm'],use_log_profile = True , z0 = 0.001) 
# Note that to be able to do a log extrapolation, the reader must have depth info

reader_roms_cnz_surface = reader_netCDF_MetOcean.Reader(data_pth + 'current_roms_cnz_surf_20050101_20050107.nc',variables_to_use = ['us','vs']) # 

reader_roms_cnz_3D = reader_netCDF_MetOcean.Reader(data_pth + 'current_roms_cnz_3D_20050101_20050107.nc',variables_to_use = ['u','v']) #

reader_swan_nzra_surface = reader_netCDF_MetOcean.Reader(data_pth + 'waves_swan_nzra_nz-nzn_20050101_20050107.nc',variables_to_use = ['hs','tp','dpm']) # 

reader_wrf_nzra1_surface = reader_netCDF_MetOcean.Reader(data_pth + 'winds_nzra1_nz_20050101_20050107.nc',variables_to_use = ['ugrd10m','vgrd10m']) # 

# Making customised landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=172.0, llcrnrlat=-42.0,
                    urcrnrlon=175, urcrnrlat=-38.0,
                    resolution='c', projection='merc') # resolution can be c (crude, the default), l (low), i (intermediate), h (high), f (full)

# In general, each variable should be input only once to avoid any confusion (expect ['x_sea_water_velocity','y_sea_water_velocity'])
# there seems to be some issues when 'dep' input twice for example
o.add_reader([reader_basemap,reader_roms_cnz_dav,reader_swan_nzra_surface,reader_wrf_nzra1_surface]) # 


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
lon = 173.7366
lat = -40.1893

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

o.seed_elements(lon, lat, radius=0, number=10,time=reader_swan_nzra_surface.start_time,
                 z=0.0, terminal_velocity = -0.001) # wind_drift_factor = 0.04 ) # age_seconds = 0,)

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
o.fallback_values['ocean_horizontal_diffusivity'] = 0.5 # specify constant ocean_horizontal_diffusivity in m2.s-1
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
o.set_config('processes:turbulentmixing', True) 
# Note that turbulentmixing include the buoyancy-related settling for now (terminal_velocity)
# if False, no vertical mixing NOR settling will occur - might need to be changed to allow pure buoyancy-drive settling and skip the vertical mixing computations?
# 
# 
o.set_config('turbulentmixing:diffusivitymodel', 'environment') # i.e. specified from model or constant
o.set_config('turbulentmixing:TSprofiles',False)
o.set_config('turbulentmixing:timestep', 1800)  
# if some ocean_vertical_diffusivity!=0, turbulentmixing:timestep should be less than 900 seconds (15min)
# if ocean_vertical_diffusivity == 0 then the vertical motion will be driven only by buoyancy set by terminal_velocity
# so this could be set the same time_step as the simulation time_step

# i.e to use a given settling velocity : terminal_velocity, with no added vertical diffusion use :
# o.fallback_values['ocean_vertical_diffusivity'] = 0.0
# The time step governing the vertical settling (and mixing if ocean_vertical_diffusivity~=0) is :
# o.set_config('turbulentmixing:timestep', 1800)  # can be set to same as "run" time step to reproduce ERcore behaviour

o.list_config()
o.list_configspec()

###############################
# RUN 
###############################
# import pdb;pdb.set_trace()
if True:
    # Running model (until end of driver data)
    o.run(time_step=1800, end_time = reader_swan_nzra_surface.start_time + timedelta(days = 7.0), outfile='opendrift_sedimentdrift3d_all_inputs.nc',time_step_output = 1800)
    # the start time is defined by seed_elements, the end_time is defined by either steps=number of step, duration = timedelta, or end_time= datetime
return

###############################
# TEST CASES = No diffusion - only one particle settling
###############################
# o.set_config('drift:scheme','euler') # or 'runge-kutta'
# o.set_config('drift:stokes_drift', False)
# o.set_config('drift:use_tabularised_stokes_drift', False)
# o.set_config('drift:tabularised_stokes_drift_fetch','25000')
# o.run(time_step=1800, end_time = reader_swan_nzra_surface.start_time + timedelta(days = 7.0), outfile='opendrift_sedimentdrift3d_euler_nostokes.nc',time_step_output = 1800)

# o.set_config('drift:scheme','runge-kutta') # or 'runge-kutta'
# o.set_config('drift:stokes_drift', False)
# o.set_config('drift:use_tabularised_stokes_drift', False)
# o.set_config('drift:tabularised_stokes_drift_fetch','25000')
# o.run(time_step=1800, end_time = reader_swan_nzra_surface.start_time + timedelta(days = 7.0), outfile='opendrift_sedimentdrift3d_runge-kutta_nostokes.nc',time_step_output = 1800)
# o.plot()

# o.set_config('drift:scheme','euler') # or 'runge-kutta'
# o.set_config('drift:stokes_drift', True)
# o.set_config('drift:use_tabularised_stokes_drift', False)
# o.set_config('drift:tabularised_stokes_drift_fetch','25000')
# o.run(time_step=1800, end_time = reader_swan_nzra_surface.start_time + timedelta(days = 7.0), outfile='opendrift_sedimentdrift3d_euler_stokes.nc',time_step_output = 1800)
# # o.plot()
# >>Stokes not compute is use_tabularised_stokes_drift:False and no [x,y] stokes components are inputs

# o.set_config('drift:scheme','euler') # or 'runge-kutta'
# o.set_config('drift:stokes_drift', True)
# o.set_config('drift:use_tabularised_stokes_drift', True)
# o.set_config('drift:tabularised_stokes_drift_fetch','25000')
# o.run(time_step=1800, end_time = reader_swan_nzra_surface.start_time + timedelta(days = 7.0), outfile='opendrift_sedimentdrift3d_euler_stokes_param.nc',time_step_output = 1800)
# o.plot()

# test windage also - can be input there > o.seed_elements(...wind_drift_factor = 0.02 ....
# o.set_config('drift:scheme','euler') # or 'runge-kutta'
# o.set_config('drift:stokes_drift', True)
# o.set_config('drift:use_tabularised_stokes_drift', True)
# o.set_config('drift:tabularised_stokes_drift_fetch','25000')
# o.run(time_step=1800, end_time = reader_swan_nzra_surface.start_time + timedelta(days = 7.0), outfile='opendrift_sedimentdrift3d_euler_stokes_param_windage4prc.nc',time_step_output = 1800)

###############################
# PLOTS / ANIMATION
###############################


# Print and plot results
# print o
# o.animation()
# o.plot()
