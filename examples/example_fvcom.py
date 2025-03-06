#!/usr/bin/env python
"""
FVCOM: Using model input from unstructured grid
===============================================
"""

from datetime import timedelta
import urllib.request as urllib_request
import numpy as np
from opendrift.readers import reader_netCDF_CF_unstructured
from opendrift.readers import reader_global_landmask
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

try:
    proj = "+proj=utm +zone=33 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    fvcom = reader_netCDF_CF_unstructured.Reader(filename = 'https://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/netcdf_unstructured_samples/AkvaplanNiva_sample_lonlat_fixed.nc', proj4 = proj)
    o.add_reader(fvcom)
    print(fvcom)
except:
    print('Thredds server not available, cannot run example')
    fvcom = None

if fvcom is not None:
    # Seed elements at defined positions, depth and time
    N = 1000
    z = -10*np.random.uniform(0, 1, N)
    o.seed_elements(lon=18.0, lat=69.8, radius=2000, number=N,
                    z=z, time=fvcom.start_time)

    #%%
    # Running model
    o.run(time_step=1800, duration=timedelta(hours=12))

    #%%
    # Print and plot results
    print(o)

    #%%
    # Animation (current as background not yet working).
    o.animation(color='z')

    o.plot(fast=True, buffer = 1.)

#%%
# .. image:: /gallery/animations/example_fvcom_0.gif

