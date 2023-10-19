#!/usr/bin/env python
"""
Drift at different depths
==========================
"""

from datetime import datetime, timedelta
import numpy as np
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Using live data from Thredds
o.add_readers_from_list([
    'https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be'])

#%%
# Adding some diffusion
o.set_config('drift:horizontal_diffusivity', 10)  # m2/s

#%%
# Seed 1000 elements at random depths
z = -np.random.rand(2000)*50
o.seed_elements(lon=4.8, lat=60.0, z=z, radius=0, number=2000,
                time=datetime.utcnow())

print(o)

#%%
# Running model
o.run(duration=timedelta(hours=24), time_step=1800)

#%%
# Plot results with lines and particles colored by depth
print(o)
o.plot(linecolor='z', buffer=.1, show_elements=False, fast=False)
o.animation(color='z', buffer=.1, fast=True)

#%%
# .. image:: /gallery/animations/example_depth_0.gif
