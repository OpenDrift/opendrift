#!/usr/bin/env python
"""
Cone seeding (Tromso)
=====================
"""

from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Using live data from Thredds
reader_arome = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')
reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader([reader_norkyst, reader_arome])

#%%
# Seed elements along cone, e.g. ship track with
# increasing uncertainty in position
latstart = 68.988911
lonstart = 16.040701
latend = 69.991446
lonend = 17.760061
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]

o.seed_elements(lon=[lonstart, lonend], lat=[latstart, latend],
                radius=[100, 500], number=10000, time=time, cone=True)

print(o)

#%%
# Adjusting some configuration
o.set_config('processes:dispersion', True)
o.set_config('processes:evaporation', False)
o.set_config('processes:emulsification', True)

#%%
# Running model for 48 hours
o.run(steps=48*2, time_step=1800, time_step_output=3600)

#%%
# Print and plot results
print(o)
o.animation()

#%%
# .. image:: /gallery/animations/example_cone_tromso_0.gif

o.plot(fast=True)
