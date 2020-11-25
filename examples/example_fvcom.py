#!/usr/bin/env python
"""
FVCOM: Using model input from unstructured grid
===============================================
"""

import numpy as np
from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.readers import reader_global_landmask
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Reader
proj = "+proj=utm +zone=33W, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
fvcom = reader_netCDF_CF_unstructured.Reader(filename = 'niva/AkvaplanNiva_sample.nc', proj4 = proj)
o.add_reader(fvcom)

# Seed elements at defined positions, depth and time
o.seed_elements(lon=17.0, lat=70.0, radius=0, number=10,
                z=np.linspace(0,-10, 10), time=fvcom.start_time)

#%%
# Running model
o.run(time_step=900)

#%%
# Print and plot results
print(o)

#%%
# Animation with current as background.
# Note that drift is also depending on wind, which is not shown.
# o.animation(background=['x_sea_water_velocity', 'y_sea_water_velocity'],
#              fast=True)

#%%
# Print and plot results
print(o)
o.plot(fast=True, buffer = 1.)
