#!/usr/bin/env python
"""
Fjord
==================================
"""

from datetime import timedelta
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.leeway import Leeway

o = Leeway(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

#%%
# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

#%%
# Making customised, full resolution landmask
reader_landmask = reader_global_landmask.Reader(
                    extent=[5.5, 6.65, 61.05, 61.21])

o.add_reader([reader_landmask, reader_norkyst, reader_arome])

#%%
# Seed elements
lat = 61.117594; lon = 6.55
#time = [reader_arome.start_time,
#        reader_arome.start_time + timedelta(hours=5)]
time = reader_arome.start_time
objType = 1  # 1: Person-in-water (PIW), unknown state (mean values)
o.seed_elements(lon, lat, radius=50, number=5000, time=time, objectType=objType)

#%%
# Running model for 66 hours, using small time step due to high resolution coastline
o.run(duration=timedelta(hours=12), time_step=300, time_step_output=3600)

#%%
# Print and plot results
print(o)
o.animation()

#%%
# .. image:: /gallery/animations/example_fjord_0.gif

o.plot()

