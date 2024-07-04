#!/usr/bin/env python
"""
SHYFEM: Using model input from unstructured grid
================================================
"""

from datetime import timedelta
import numpy as np
from opendrift.readers.unstructured import shyfem
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

o.set_config('general:coastline_action', 'previous')

shyfem = shyfem.Reader('https://iws.ismar.cnr.it/thredds/dodsC/emerge/shyfem_unstructured_adriatic.nc')
o.add_reader(shyfem)
print(shyfem)

# Seed elements at defined positions, depth and time
N = 1000
z = -15*np.random.uniform(0, 1, N)
o.seed_elements(lon=12.4, lat=45.25, radius=1000, number=N,
                z=z, time=shyfem.start_time)

#%%
# Running model
o.run(time_step=1800, duration=timedelta(hours=12))

#%%
# Print and plot results
print(o)

#%%
# Animations
o.animation(color='z', markersize=5, corners=[12.2, 12.6, 45.1, 45.5])

o.animation_profile(color='z', markersize=5)


#%%
# .. image:: /gallery/animations/example_shyfem_0.gif
# .. image:: /gallery/animations/example_shyfem_1.gif
o.plot(fast=True, buffer = 1., corners=[12.2, 12.6, 45.1, 45.5])
