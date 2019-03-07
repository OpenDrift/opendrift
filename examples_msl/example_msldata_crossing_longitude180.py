#!usrbinenv python

#  This script is based on httpsgithub.comOpenDriftopendrifttreemasterexamplesexample_grid.py
#  The aim is to test the use of reader_netCDF_MetOcean which allows using MetOcean data in OpenDrift
#   
#  This is also used as a reference case to compare Opendrift trajectories with ERcore trajectories
# 

import numpy as np
from datetime import datetime, timedelta
from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_netCDF_MetOcean
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers.interpolation import ReaderBlock

###############################
# MODEL SELECTION
###############################

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
###############################
# READERS
###############################
reader_cfsr = reader_netCDF_MetOcean.Reader('E:/metocean/0440_Tawhaki1_drillcuttings/opendrift_modelling/cfsr_ocean_nz_20040101.nc')
reader_roms_nz = reader_netCDF_MetOcean.Reader('E:/metocean/0440_Tawhaki1_drillcuttings/opendrift_modelling/cfsr_ocean_nz_20040101.nc')#BL=[165.03 -47.97] UR=[179.96 -33.02]

# Making customised landmask (Basemap)

reader_basemap = reader_basemap_landmask.Reader(
                llcrnrlon=160.25, llcrnrlat=-50.75, 
                urcrnrlon=189.75, urcrnrlat=-33.02,   
                resolution='l', projection='merc') # resolution can be c (crude, the default), l (low), i (intermediate), h (high), f (full)    

# o.add_reader([reader_basemap, reader_norkyst])
o.add_reader([reader_basemap,reader_cfsr])

o.fallback_values['x_wind'] = 10 # strong westerly wind to make particle drift past 180degT and check behaviour
o.fallback_values['y_wind'] = 0
# o.fallback_values['x_wind'] = 0 # strong westerly wind to make particle drift past 180degT and check behaviour
# o.fallback_values['y_wind'] = 0

###############################
# PARTICLE SEEDING
###############################

# Seeding some particles
# lons = np.linspace(3.5, 5.0, 100)
# lats = np.linspace(60, 61, 100)
# 
#  Point release
lon = 179.9; lat = -48; 
#  Rectangle release
lons = np.linspace(179.0,181.0, 10) 
lats = np.linspace(-45.5,-45.5, 10)
# these lon/lat spread across both roms_nz and cfsr data
depth = np.linspace(0.,0., 100) # surface=0, negative down

lons, lats = np.meshgrid(lons, lats)

lons = lons.ravel()
lats = lats.ravel()

#cfsr
v = reader_cfsr.get_variables(['y_sea_water_velocity', 'x_sea_water_velocity'], time=reader_cfsr.start_time, x=lons, y=lats, z=depth, block=True)
b = ReaderBlock(v, interpolation_horizontal='linearND')
env, prof = b.interpolate(lons,lats, depth, ['y_sea_water_velocity', 'x_sea_water_velocity'])
print env
# roms_nz
v1 = reader_roms_nz.get_variables(['y_sea_water_velocity', 'x_sea_water_velocity'], time=reader_cfsr.start_time, x=lons, y=lats, z=depth, block=True)
b1 = ReaderBlock(v1, interpolation_horizontal='linearND')
env1, prof1 = b1.interpolate(lons,lats, depth, ['y_sea_water_velocity', 'x_sea_water_velocity'])
print env1
#
import pdb;pdb.set_trace()
# Seed oil elements at defined position and time

o.seed_elements(lons, lats, radius=0, number=100,
                time = reader_cfsr.start_time, wind_drift_factor = 0.5) # large wind_drift_factor to make particles drift past 180 degT

###############################
# PHYSICS
###############################
# o.set_config('driftcurrent_uncertainty', .1)

###############################
# RUN 
###############################
# Running model (until end of driver data)
o.run(time_step=1800, end_time = reader_cfsr.start_time + timedelta(days = 1), outfile='opendrift_test_longitude180_crossing.nc',time_step_output = 1800)  
# the start time is define by seed_elements, the end_time is defined by either steps=number of step, duration = timedelta, or end_time= datetime

###############################
# PLOTS  ANIMATION
###############################
# Print and plot results
print o
o.animation()
o.plot()
