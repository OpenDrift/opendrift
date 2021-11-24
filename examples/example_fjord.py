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
# Add readers for wind and current
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
o.add_reader([reader_norkyst, reader_arome])

#%%
# Seed elements
#time = [reader_arome.start_time,
#        reader_arome.start_time + timedelta(hours=5)]
time = reader_arome.start_time
object_type = 1  # 1: Person-in-water (PIW), unknown state (mean values)
o.seed_elements(lon=6.55, lat=61.117594, radius=50, number=5000, time=time, object_type=object_type)

#%%
# Running model for 12 hours, using small time step due to high resolution coastline
o.run(duration=timedelta(hours=12), time_step=300, time_step_output=3600)

#%%
# Print and plot results
print(o)
o.animation()

#%%
# .. image:: /gallery/animations/example_fjord_0.gif

o.plot()

