#!/usr/bin/env python
"""
Cone seeding
=====================
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Using live data from Thredds
o.add_readers_from_list([
    'https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be'])

#%%
# Seed elements along cone, e.g. ship track with
# increasing uncertainty in position
latstart = 68.988911
lonstart = 16.040701
latend = 69.991446
lonend = 17.760061
time = [datetime.utcnow(), datetime.utcnow() + timedelta(hours=12)]

o.seed_cone(lon=[lonstart, lonend], lat=[latstart, latend],
            oiltype='EKOFISK', radius=[100, 800], number=10000, time=time)

print(o)

#%%
# Adjusting some configuration
o.set_config('processes:dispersion', True)
o.set_config('processes:evaporation', False)
o.set_config('processes:emulsification', True)
o.disable_vertical_motion()
#o.set_config('drift:vertical_mixing', False)
#o.set_config('drift:vertical_mixing', False)

#%%
# Running model for 24 hours
o.run(steps=24*2, time_step=1800, time_step_output=3600)

#%%
# Print and plot results
print(o)
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_cone_0.gif

o.plot(fast=True)
