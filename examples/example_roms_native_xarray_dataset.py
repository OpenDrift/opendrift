#!/usr/bin/env python
"""
ROMS native reader
==================================
"""

import numpy as np
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift
import xarray as xr

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

#%%
filename = (o.test_data_folder() +
    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

# Creating and adding reader for Nordic 4km current dataset

ds = xr.open_dataset(filename, decode_times=False)
nordic_native = reader_ROMS_native.Reader(ds)
o.add_reader(nordic_native)

#%%
# Seed elements at defined positions, depth and time
o.seed_elements(lon=12.0, lat=68.3, radius=0, number=10,
                z=np.linspace(0, -150, 10), time=nordic_native.start_time)

#%%
# Running model
o.run(time_step=3600)

#%%
# Print and plot results, with lines colored by particle depth
print(o)
o.plot(linecolor='z', fast=True)
#o.animation()
