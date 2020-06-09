#!/usr/bin/env python
"""
Cod egg
=============
"""

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.pelagicegg import PelagicEggDrift
import datetime

o = PelagicEggDrift(loglevel=20)  # Set loglevel to 0 for debug information

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norshelf = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea_norshelf_his_ZDEPTHS_agg')

o.add_reader([reader_norkyst, reader_norshelf])

#%%
# spawn NEA cod eggs at defined position and time
time = datetime.datetime.today() - datetime.timedelta(1)
o.seed_elements(14. , 68.1, z=-40, radius=2000, number=500,
                time=time, diameter=0.0014, neutral_buoyancy_salinity=31.25)
o.seed_elements(12.5, 68., z=-40, radius=2000, number=500,
                time=time, diameter=0.0014, neutral_buoyancy_salinity=31.25)
o.seed_elements(13.5, 68., z=-40, radius=2000, number=500,
                time=time, diameter=0.0014, neutral_buoyancy_salinity=31.25)
o.seed_elements(13., 67.8, z=-40, radius=2000, number=500,
                time=time, diameter=0.0014, neutral_buoyancy_salinity=31.25)



#%%
# Adjusting some configuration
o.set_config('processes:turbulentmixing', True)
#o.set_config('turbulentmixing:diffusivitymodel', 'windspeed_Sundby1983') # windspeed parameterization for eddy diffusivity
o.set_config('turbulentmixing:diffusivitymodel', 'environment') # use eddy diffusivity from ocean model
#%%
# Vertical mixing requires fast time step
o.set_config('turbulentmixing:timestep', 60.) # seconds

#%%
# Running model (until end of driver data)
o.run(steps=96, time_step=3600)

#%%
# Print and plot results
print(o)

o.plot(fast=True)
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_codegg_0.gif

o.plot_vertical_distribution()
