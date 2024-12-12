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
# Adjusting some configuration
o.set_config('processes:dispersion', True)
o.set_config('processes:evaporation', False)
o.set_config('processes:emulsification', True)

#%%
# Seed elements along cone, e.g. ship track with
# increasing uncertainty in position
latstart = 68.988911
lonstart = 16.040701
latend = 69.991446
lonend = 17.760061
time = [datetime.utcnow(), datetime.utcnow() + timedelta(hours=12)]
o.seed_cone(lon=[lonstart, lonend], lat=[latstart, latend],
            oil_type='EKOFISK', radius=[100, 800], number=10000, time=[time])

print(o)

#%%
# Running model for 24 hours
o.run(steps=24*2, time_step=1800, time_step_output=3600)

#%%
# Print and plot results
print(o)

#%%
# Add text label on the map
text = [{'s': 'Senja', 'x': 17.3, 'y': 69.3, 'fontsize': 20, 'color': 'g',
         'backgroundcolor': 'white', 'bbox': dict(facecolor='white', alpha=0.8), 'zorder': 1000}]

o.animation(fast=False, ocean_color='skyblue', land_color='burlywood', text=text)

#%%
# .. image:: /gallery/animations/example_cone_0.gif

o.plot(fast=True, ocean_color='skyblue', land_color='dimgray', text=text)
