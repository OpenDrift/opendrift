#!/usr/bin/env python

#  This script is based on https://github.com/OpenDrift/opendrift/tree/master/examples_msl/example_grid_msldata_settling.py
#  The aim is to show how to:
#  - specify several readers, with (u,v) pairs, adding their velocities for the simulations
#  - set of primary/secondary readers to be used depending on particle positions
#   
#  This is also used as a reference case to compare Opendrift trajectories with ERcore trajectories
# 
#
#  This example show how to add readers for the currents
 # for example the tidal and residual components
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
# depth
# reader_roms_bob_depth = reader_netCDF_MetOcean.Reader('C:/metocean/bob_cons.nc',variables_to_use = ['dep']) # i.e. only use the 'dep' variable in that file
# reader_roms_bob_depth.always_valid = True # this is useful if the time vector associated with the bathy does not fit with the simulation time 
#use always_valid for variables that can be used at all times, such as a bathy 

# residual currents
reader_roms_res3D = reader_netCDF_MetOcean.Reader('F:/metocean/0000_OMV_drillcuttings/opendrift_modelling/flow_fields/roms_cnz-residual_3D_20040101.nc',variables_to_use = ['uo','vo']) # 
# contents of reader can be check by calling reader_roms_bob_res3D - this  will show which variables it includes, timing, levels, extents etc.. etc..
# tidal currents
reader_roms_tide = reader_netCDF_MetOcean.Reader('F:/metocean/0000_OMV_drillcuttings/opendrift_modelling/flow_fields/nz_tide_20040101.nc',variables_to_use = ['ut','vt']) # 

reader_cur_total = reader_roms_tide + reader_roms_res3D

reader_basemap = reader_basemap_landmask.Reader(
                llcrnrlon=172.65, llcrnrlat=-39.6262, 
                urcrnrlon=174.65, urcrnrlat=-37.6262,   
                resolution='c', projection='merc') # resolution can be c (crude, the default), l (low), i (intermediate), h (high), f (full)


o.add_reader([reader_basemap,reader_cur_total]) #
# note the order in which 'readers' are specified matters - the entered first will be used as first choice, second one as fallback etc... 

# import pdb;pdb.set_trace()
o.list_environment_variables()

###############################
# PARTICLE SEEDING
###############################

#  Point release
lon = 173.652299949710
lat = -38.626275006462530

# Seed oil elements at defined position and time

o.seed_elements(lon, lat, radius=0, number=1000,time=reader_roms_res3D.start_time,
                 z=0.0, terminal_velocity = -0.001) #, wind_drift_factor = 0, age_seconds = 0,)

# specific element variable such as terminal_velocity, can be specified here. 
# terminal_velocity>0 particle moves up, terminal_velocity<0 particle moves down
# 
# General variables to be defined lon,lat,z,time ..look in each /elements/element*.py for more properties
# 
# 
# import pdb;pdb.set_trace()

###############################
# PHYSICS
###############################
# Default

###############################
# RUN CASE WITH ADDED CURRENTS : TIDE+RES
###############################

# Running model (until end of driver data)
o.run(time_step=1800, end_time = reader_roms_res3D.start_time+timedelta(days=3.0), outfile='opendrift_tide_plus_res.nc',time_step_output = 1800)
# the start time is defined by seed_elements, the end_time is defined by either steps=number of step, duration = timedelta, or end_time= datetime
# o.plot()

#############################################
# RUN CASES WITH PURE TIDE AND PURE RESIDUALS
#############################################
# 
# TIDE
# 
o1 = SedimentDrift3D(loglevel=0)  # Set loglevel to 0 for debug information
del reader_roms_tide
reader_roms_tide = reader_netCDF_MetOcean.Reader('F:/metocean/0000_OMV_drillcuttings/opendrift_modelling/flow_fields/nz_tide_20040101.nc',variables_to_use = ['ut','vt']) # 
# reader_roms_tide > MUST BE RE-INITIALZED otherwise it will keep adding the residuals...
o1.add_reader([reader_basemap,reader_roms_tide]) #

o1.seed_elements(lon, lat, radius=0, number=1000,time=reader_roms_res3D.start_time,
                 z=0.0, terminal_velocity = -0.001) #, wind_drift_factor = 0, age_seconds = 0,)

o1.run(time_step=1800, end_time = reader_roms_res3D.start_time+timedelta(days=3.0), outfile='opendrift_tide_only.nc',time_step_output = 1800)
# 
# RESIDUAL
# 
o2 = SedimentDrift3D(loglevel=0)  # Set loglevel to 0 for debug information
o2.add_reader([reader_basemap,reader_roms_res3D]) #

o2.seed_elements(lon, lat, radius=0, number=1000,time=reader_roms_res3D.start_time,
                 z=0.0, terminal_velocity = -0.001) #, wind_drift_factor = 0, age_seconds = 0,)

o2.run(time_step=1800, end_time = reader_roms_res3D.start_time+timedelta(days=3.0), outfile='opendrift_res_only.nc',time_step_output = 1800)


###############################
# PLOTS / ANIMATION
###############################
import pdb;pdb.set_trace()
# Print and plot results
# o.animation()
o.plot()
o1.plot()
o2.plot()
o.animation(compare=o1,legend=['tide+res', 'tide only'])
o.animation(compare=o2, legend=['tide+res', 'residual only'])