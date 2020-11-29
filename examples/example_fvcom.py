#!/usr/bin/env python
"""
FVCOM: Using model input from unstructured grid
===============================================
"""

from datetime import timedelta
import numpy as np
from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.readers import reader_global_landmask
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Reader
proj = "+proj=utm +zone=33W, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
fvcom = reader_netCDF_CF_unstructured.Reader(filename = 'https://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/netcdf_unstructured_samples/AkvaplanNiva_sample_lonlat_fixed.nc', proj4 = proj)
o.add_reader(fvcom)
print(fvcom)

# Seed elements at defined positions, depth and time
N = 1000
z = -10*np.random.uniform(0, 1, N)
o.seed_elements(lon=17.0, lat=70.0, radius=5000, number=N,
                z=z, time=fvcom.start_time)

#%%
# Running model
o.run(time_step=1800, duration=timedelta(hours=12), outfile='fvcom.nc')

#%%
# Print and plot results
print(o)

#%%
# Animation (current as background not yet working).
# Note that drift is also depending on wind, which is not shown.
o.animation(filename='fvcom.mp4', color='z')
#o.animation(background=['x_sea_water_velocity', 'y_sea_water_velocity'],
#            fast=True)

#%%
# Print and plot results
print(o)
o.plot(fast=True, buffer = 1.)
