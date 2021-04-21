#!/usr/bin/env python
"""
SHYFEM: Using model input from unstructured grid
===============================================
"""

from datetime import timedelta
import numpy as np
from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.readers.unstructured import shyfem
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

shyfem = shyfem.Reader('https://iws.ismar.cnr.it/thredds/dodsC/emerge/shyfem_unstrct_ADRIA_20210408.nc')
o.add_reader(shyfem)
print(shyfem)

# Seed elements at defined positions, depth and time
N = 1000
z = -10*np.random.uniform(0, 1, N)
o.seed_elements(lon=13., lat=40., radius=2000, number=N,
                z=z, time=shyfem.start_time)

#%%
# Running model
o.run(time_step=1800, duration=timedelta(hours=12))

#%%
# Print and plot results
print(o)

#%%
# Animation (current as background not yet working).
o.animation(color='z')

#%%
# .. image:: /gallery/animations/example_fvcom_0.gif

o.plot(fast=True, buffer = 1.)
