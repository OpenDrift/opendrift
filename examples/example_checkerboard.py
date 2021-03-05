#!/usr/bin/env python
"""
Checkerboard
============
"""

import numpy as np
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([reader_norkyst])

o.set_config('drift:vertical_mixing', False)
#%%
# Seeding particles in a checkerboard pattern
di = 5 # Horizontal number of particles per square
dj = 5 # Vertical number of particles per square
lons = np.linspace(3.5, 5.0, 100)
lats = np.linspace(60, 61, 100)

ii = np.arange(len(lons))//di
jj = np.arange(len(lats))//dj
ii, jj = np.meshgrid(ii, jj)
board = (ii+jj)%2 > 0

lons, lats = np.meshgrid(lons, lats)
lons = lons[board].ravel()
lats = lats[board].ravel()

o.seed_elements(lons, lats, radius=0, number=5000,
                time=reader_norkyst.start_time)

#%%
# Running model
o.run(steps=66*2, time_step=1800)

#%%
# Print and plot results
print(o)
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_checkerboard_0.gif

o.plot()
