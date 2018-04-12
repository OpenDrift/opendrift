#!/usr/bin/env python

# Simple test script showing the use of the SedimentDrift3D module. 
# The script is based on example_openoil3d.py

from datetime import datetime, timedelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.sedimentdrift3D import SedimentDrift3D

o = SedimentDrift3D(loglevel=0) # SedimentDrift3D is based on OpenDrift3DSimulation base class

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([reader_norkyst, reader_arome])
# NOTE : in most case, we will need to have the variable 'sea_floor_depth_below_sea_level'
# defined by one of the reader in order to evaluate if particles have touched the seabed or not.
# This is not the case here as we are simply re-using some existing test data, so particle will 
# settle indefinitely.


# Seeding some particles
lon = 4.9; lat = 60.1; # Outside Bergen
time = reader_arome.start_time

# Seed sediment elements at defined position and time
# the sediment settling velocity is defined by terminal_velocity [m/s]
# terminal_velocity<0 : negatively buoyant, i.e. settle towards the seabed
# terminal_velocity>0 : positively buoyant, i.e. rise towards the surface
o.seed_elements(lon, lat, radius=10, number=1000, time=time, z=0,terminal_velocity = -0.001)

# Adjusting some configuration
# ocean_horizontal_diffusivity in m2.s-1 - new feature implemented in SedimentDrift3D, could be extended to basemodel.py
o.fallback_values['ocean_horizontal_diffusivity'] = 0.1
# specify constant ocean_vertical_diffusivity in m2.s-1 
o.fallback_values['ocean_vertical_diffusivity'] = 0.0001 

#processes
o.set_config('processes:verticaladvection' , False) # no vertical current available, so no vertical advection
# New process added to SedimentDrift3D : resuspension
# This still needs to be implemented/tested but will eventually allow checking bed shear stress at particle position 
# and evaualte if resuspension is possible based on a "critical" shear stress
# already False be default but mentioned here for reference
o.set_config('processes:resuspension',False) 
o.set_config('processes:turbulentmixing', True) 
o.set_config('turbulentmixing:timestep',  5)

# Running model (until end of driver data)
o.run(steps=4*40, time_step=900,
      time_step_output=3600,
      outfile='sediment3d.nc')

# Print and plot results
print o
o.animation()
o.plot()
