#!/usr/bin/env python
"""
Sediment simulation (with settling on seafloor)

This script is re-using the same data/release details as example_long_oilspill_seafloor.py

==================================
"""

from datetime import timedelta

from opendrift.readers import reader_netCDF_CF_generic
# from opendrift.models.openoil3D import OpenOil3D
from opendrift.models.sedimentdrift3D import SedimentDrift3D # import sedimentDrift3D module
import numpy as np

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

###############################
# CONFIG
###############################
o.set_config('environment:fallback:x_wind', 0)
o.set_config('environment:fallback:y_wind', 0)
# Adjusting some configuration
# diffusion -
Kxy = 0.1 # m2/s-1
Kz = 0.0001 # m2/s-1
o.set_config('drift:horizontal_diffusivity', Kxy)     # Difference from first run
# specify constant ocean_vertical_diffusivity in m2.s-1 
o.set_config('environment:fallback:ocean_vertical_diffusivity',Kz)
# processes
o.set_config('drift:vertical_advection' , False) # no vertical current available, so no vertical advection
o.set_config('drift:vertical_mixing', True) 
o.set_config('vertical_mixing:timestep',  5)

# New processes added to SedimentDrift3D : resuspension, and online computation of terminal_velocity
o.set_config('drift:resuspension',False) 
o.set_config('processes:update_terminal_velocity',True) 
o.set_config('processes:terminal_velocity_method','van_rijn1984')  # ['dietrich','stokes','zhiyao','ahrens','komar','van_rijn1984','soulsby1997']

# Seeding some particles
lon = 4.5; lat = 62.0
time = [reader_norkyst.start_time,
        reader_norkyst.end_time]

# Seed sediment elements at defined position and time
# the sediment settling velocity is defined by terminal_velocity [m/s]
# terminal_velocity<0 : negatively buoyant, i.e. settle towards the seabed
# terminal_velocity>0 : positively buoyant, i.e. rise towards the surface

if False:
    # constant terminal_velocity input 
    o.set_config('processes:update_terminal_velocity',False)       
    o.seed_elements(lon = lon, 
            lat = lat, 
            z=-1.0, 
            radius=5, 
            number=3000, 
            time=time,
            terminal_velocity = -0.03)

if True :
    o.set_config('processes:update_terminal_velocity',True) 
    o.set_config('processes:terminal_velocity_method','van_rijn1984')  # ['dietrich','stokes','zhiyao','ahrens','komar','van_rijn1984','soulsby1997']
    d50_array = 125e-6 * np.ones(3000) + 100e-6 * np.random.rand(3000) # mimic a PSD, centered around 125 microns
    o.seed_elements(lon = lon, 
                    lat = lat, 
                    z=-1.0, 
                    radius=5, 
                    number=3000, 
                    d50 = d50_array, # median grain size - could be an array with PSD, with some randomness
                    terminal_velocity = -0.001,
                    critical_shear_stress = 0.01,
                    time=time)

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

