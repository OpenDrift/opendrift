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
# from opendrift.models.sedimentdrift3D import SedimentDrift3D
# from opendrift.models.oceandrift3D import OceanDrift3D
from opendrift.models.openoil3D import OpenOil3D

###############################
# MODEL SELECTION
###############################
o = OpenOil3D(loglevel=0, weathering_model = 'default') #weathering_model = 'default' # Using Oil information from \opendrift\models\oilprop.dat

print o.oiltypes  #print available oil types

o.max_speed = 5.0 # to make sure that netcdf frames used for interpolation are large enough

###############################
# READERS
###############################
data_pth = 'C:/Users/simon/Google Drive/R&D/OpenDrift/data_testing/'
# DAV currents
# reader_roms_cnz_dav = reader_netCDF_MetOcean.Reader(data_pth + 'current_roms_cnz_surf_20050101_20050107_new.nc',variables_to_use = ['dep','um','vm'],use_log_profile = True , z0 = 0.001) 
# Note that to be able to do a log extrapolation, the reader must have depth info

# SURFACE currents
reader_roms_cnz_surface = reader_netCDF_MetOcean.Reader(data_pth + 'current_roms_cnz_surf_20050101_20050107_new.nc',variables_to_use = ['dep','us','vs']) # 
# 3D currents
# reader_roms_cnz_3D = reader_netCDF_MetOcean.Reader(data_pth + 'current_roms_cnz_3D_20050101_20050107_new.nc',variables_to_use = ['temp']) #
# WATER TEMP - needs to be 3D to work
reader_roms_cnz_sea_temp = reader_netCDF_MetOcean.Reader(data_pth + 'current_roms_cnz_3D_20050101_20050107_new.nc',variables_to_use = ['temp','salt']) #
# WAVES
reader_swan_nzra_surface = reader_netCDF_MetOcean.Reader(data_pth + 'waves_swan_nzra_nz-nzn_20050101_20050107.nc',variables_to_use = ['hs','tp','dpm']) # 
# WINDS
reader_wrf_nzra1_surface = reader_netCDF_MetOcean.Reader(data_pth + 'winds_nzra1_nz_20050101_20050107.nc',variables_to_use = ['ugrd10m','vgrd10m']) # 

# Making customised landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=172.0, llcrnrlat=-42.0,
                    urcrnrlon=175, urcrnrlat=-38.0,
                    resolution='c', projection='merc') # resolution can be c (crude, the default), l (low), i (intermediate), h (high), f (full)

# In general, each variable should be input only once to avoid any confusion (expect ['x_sea_water_velocity','y_sea_water_velocity'])
# there seems to be some issues when 'dep' input twice for example
o.add_reader([reader_basemap,reader_roms_cnz_surface,reader_roms_cnz_sea_temp,reader_swan_nzra_surface,reader_wrf_nzra1_surface]) # 

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

# Seed oil elements at defined position and time

# o.seed_elements(lons, lats, radius=0, number=10000,
#                 time=reader_roms_cnz.start_time)

o.seed_elements(lon, lat, 
                radius=10,
                number=100,
                time=reader_roms_cnz_surface.start_time,
                z=0.0,
                oiltype='MARIN DIESEL', # must be picked from available o.oil_types
                # m3_per_hour =***) # 
                # diameter = *** , to force diameter of entrained droplets
                wind_drift_factor = 0.04 ) # all possible specs are in openoil.py

###############################
# PHYSICS
###############################

# these will list all possible options for that model
# diffusion - constant in that example
o.fallback_values['ocean_horizontal_diffusivity'] = 0.0 # specify constant ocean_horizontal_diffusivity in m2.s-1
o.fallback_values['ocean_vertical_diffusivity'] = 0.0000 # specify constant ocean_vertical_diffusivity in m2.s-1
o.fallback_values['sea_water_temperature'] = 15.0 # specify constant sea_water_temperature in degC
o.fallback_values['sea_water_salinity'] = 35.0 # specify constant sea_water_salinity in ppm

# drift
o.set_config('drift:scheme','runge-kutta') # or 'runge-kutta'
o.set_config('drift:current_uncertainty', 0.05)
o.set_config('drift:wind_uncertainty', 1.0)
o.set_config('drift:stokes_drift', False) # not available in wave files
o.set_config('drift:use_tabularised_stokes_drift', False) # not using estimations from wind for now.
o.set_config('drift:tabularised_stokes_drift_fetch','25000')

o.set_config('processes:evaporation',  True)
o.set_config('processes:emulsification',  True)
o.set_config('processes:turbulentmixing',  True)
o.set_config('turbulentmixing:diffusivitymodel', 'environment') # i.e. specified from model or constant
o.set_config('turbulentmixing:TSprofiles',False)
o.set_config('turbulentmixing:timestep', 60)  
#processes
o.set_config('processes:turbulentmixing', True) 
o.set_config('processes:verticaladvection' , False) # no vertical current available, so no vertical advection
# o.set_config('processes:resuspension',False) # already False be default but just for reference 

o.list_config()
o.list_configspec()

###############################
# RUN 
###############################
# Running model
o.run(time_step = 1800, end_time = reader_roms_cnz_surface.start_time + timedelta(days = 7.0), outfile='opendrift_oil3d_all_inputs.nc',time_step_output = 1800)
# the start time is defined by seed_elements, the end_time is defined by either steps=number of step, duration = timedelta, or end_time= datetime

###############################
# PLOTS / ANIMATION
###############################

# Print and plot results
print o
o.animation()
o.plot()
