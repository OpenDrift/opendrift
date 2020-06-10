#!/usr/bin/env python
"""
Back and forth
==============
"""

import os
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift


ncfile = 'backandforth.nc'

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader(reader_norkyst)

#%%
# Forward run
# Seeding some particles
lon = 4.2; lat = 60.1;
time = reader_norkyst.start_time
o.seed_elements(lon, lat, radius=1000, number=100, time=time)

o.run(steps=50*4, time_step=900, outfile=ncfile)

#%%
# Print and plot results
print(o)
o.plot(buffer=.2, Fast=True)

##%%
## Backward run
o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Import forward run, and seed elements at final positions and time
o.io_import_file(ncfile)
elements_final = o.elements
time_final = o.time
del o
o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o.add_reader(reader_norkyst)
o.schedule_elements(elements_final, time_final)

#%%
# Running model backwards from end of forward simulation
o.run(steps=50*4, time_step=-900)

#%%
# Print and plot results
print(o)
o.plot(buffer=.2, fast=True)
os.remove(ncfile)

##%%
# Compare plots forward and backward
