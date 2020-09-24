#!/usr/bin/env python
"""
Drift at different depths
==========================
"""

from datetime import datetime, timedelta
import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Using live data from Thredds
o.add_readers_from_list([
    'https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be'])

#%%
# Seed 1000 elements at random depths
z = -np.random.rand(1000)*50
o.seed_elements(lon=4.8, lat=60.0, z=z, radius=0, number=1000,
                time=datetime.utcnow())

print(o)

#%%
# Adding some diffusion
o.set_config('drift:current_uncertainty', .1)
o.set_config('drift:wind_uncertainty', 2)

#%%
# Running model
o.run(duration=timedelta(hours=24), time_step=1800, outfile='openoil.nc',
      export_buffer_length=5)  # Writing to netCDF file every 5 time steps

#%%
# Plot results with lines and particles colored by depth
print(o)
o.plot(linecolor='z', buffer=.1, show_particles=False, fast=True)
o.animation(color='z', buffer=.1, fast=True)

#%%
# .. image:: /gallery/animations/example_depth_0.gif
