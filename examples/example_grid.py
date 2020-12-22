#!/usr/bin/env python
"""
Grid
=============
"""

import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader(reader_norkyst)
#o.set_config('environment:fallback:land_binary_mask', 0)
o.set_config('drift:vertical_mixing', False)

#%%
# Seeding some particles
lons = np.linspace(3.5, 5.0, 100)
lats = np.linspace(60, 61, 100)
lons, lats = np.meshgrid(lons, lats)
lons = lons.ravel()
lats = lats.ravel()
o.seed_elements(lons, lats, radius=0, number=10000,
                time=reader_norkyst.start_time)

#%%
# Running model
o.run(steps=60*2, time_step=1800, time_step_output=1800)

#%%
# Print and plot results
print(o)
o.animation(fast=False, corners=[3.5, 5.5, 59.9, 61.2])

#%%
# .. image:: /gallery/animations/example_grid_0.gif
