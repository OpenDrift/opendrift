#!/usr/bin/env python
"""
Same as example_sedimentdrif3d.py, but assuming three sediment class with different settling 
velocities

Sediment simulation (with settling on seafloor)

This script is re-using the same data/release details as example_long_oilspill_seafloor.py

==================================
"""

from datetime import timedelta

from opendrift.readers import reader_netCDF_CF_generic
# from opendrift.models.openoil3D import OpenOil3D
from opendrift.models.sedimentdrift3D import SedimentDrift3D # import sedimentDrift3D module


o = SedimentDrift3D(loglevel=0)  # Set loglevel to 0 for debug information

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

# The netcdf file has data for the following variables : 
#   sea_floor_depth_below_sea_level
#   sea_water_salinity
#   sea_water_temperature
#   x_sea_water_velocity
#   y_sea_water_velocity
# 
# 	NOTE : we need to have the variable 'sea_floor_depth_below_sea_level'
# 	defined by one of the reader in order to evaluate if particles have touched the seabed or not.

o.add_reader([reader_norkyst])

o.set_config('environment:fallback:x_wind', 0)
o.set_config('environment:fallback:y_wind', 0)
# Adjusting some configuration
# diffusion -
Kxy = 0.1 # m2/s-1
Kz = 0.0001 # m2/s-1

o.set_config('drift:horizontal_diffusivity', Kxy)     # Difference from first run
# specify constant ocean_vertical_diffusivity in m2.s-1 
o.set_config('environment:fallback:ocean_vertical_diffusivity',Kz)

# Seeding some particles
lon = 4.5; lat = 62.0
time = [reader_norkyst.start_time,
        reader_norkyst.end_time]

# Seed sediment elements at defined position and time
# the sediment settling velocity is defined by terminal_velocity [m/s]
# terminal_velocity<0 : negatively buoyant, i.e. settle towards the seabed
# terminal_velocity>0 : positively buoyant, i.e. rise towards the surface
o.seed_elements(lon, lat, z=-1.0, radius=5, number=3000, time=time,terminal_velocity = -0.03)
o.seed_elements(lon, lat, z=-1.0, radius=5, number=3000, time=time,terminal_velocity = -0.003)
o.seed_elements(lon, lat, z=-1.0, radius=5, number=3000, time=time,terminal_velocity = -0.0003)

# processes
o.set_config('drift:vertical_advection' , False) # no vertical current available, so no vertical advection
# New process added to SedimentDrift3D : resuspension
# This still needs to be implemented/tested but will eventually allow checking bed shear stress at particle position 
# and evaualte if resuspension is possible based on a "critical" shear stress
# already False be default but mentioned here for reference
o.set_config('drift:resuspension',False) 
o.set_config('drift:vertical_mixing', True) 
o.set_config('vertical_mixing:timestep',  5)

# Running model
# o.run(steps=5*12*2, time_step=300, time_step_output=300)
o.run(end_time=reader_norkyst.end_time, time_step=300, time_step_output=300)

# Print and plot results
print(o)
o.plot(fast=True)
o.animation_profile() #filename='sediment_settling_profile.gif')
o.animation(color='z', background='sea_floor_depth_below_sea_level')#,filename='sediment_settling_2D.gif',)
#%%
# .. image:: /gallery/animations/sediment_settling_from_surface.gif

