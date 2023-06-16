#!/usr/bin/env python
"""
Grid time
=============
"""

from datetime import timedelta
import numpy as np
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

# Seeding at a grid at regular interval
o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

#%%
# Landmask
reader_landmask = reader_global_landmask.Reader()

o.add_reader([reader_landmask, reader_norkyst])

#%%
# Seeding some particles
lons = np.linspace(4.4, 4.6, 10)
lats = np.linspace(60.0, 60.1, 10)
lons, lats = np.meshgrid(lons, lats)
lons = lons.ravel()
lats = lats.ravel()

#%%
# Seed oil elements on a grid at regular time interval
start_time = reader_norkyst.start_time
time_step = timedelta(hours=6)
num_steps = 10
for i in range(num_steps+1):
    o.seed_elements(lons, lats, radius=0, number=100,
                    time=start_time + i*time_step)

#%%
# Running model for 60 hours
o.run(steps=60*4, time_step=900, time_step_output=3600)

#%%
# Print and plot results
print(o)
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_grid_time_0.gif
