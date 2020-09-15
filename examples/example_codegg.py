#!/usr/bin/env python
"""
Cod egg
=============
"""

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.pelagicegg import PelagicEggDrift
from datetime import datetime

o = PelagicEggDrift(loglevel=20)  # Set loglevel to 0 for debug information

# Forcing with Topaz ocean model and MEPS atmospheric model
o.add_readers_from_list([
    'https://thredds.met.no/thredds/dodsC/topaz/dataset-topaz4-arc-unmasked-be',
    'https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc'])

#%%
# spawn NEA cod eggs at defined position and time
time = datetime.utcnow()
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
o.set_config('drift:vertical_mixing', True)
#o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983') # windspeed parameterization for eddy diffusivity
o.set_config('vertical_mixing:diffusivitymodel', 'environment') # use eddy diffusivity from ocean model
#%%
# Vertical mixing requires fast time step
o.set_config('vertical_mixing:timestep', 60.) # seconds

#%%
# Running model
o.run(steps=96, time_step=3600)

#%%
# Print and plot results.
# At the end the wind vanishes, and eggs come to surface
print(o)

o.plot(fast=True)
o.animation(fast=True, color='z')

#%%
# .. image:: /gallery/animations/example_codegg_0.gif

#%% Interactive slider (not working in browser)
o.plot_vertical_distribution()
