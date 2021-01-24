#!/usr/bin/env python
"""
Fish Eggs and Larvae
====================
"""

from datetime import datetime, timedelta
from opendrift.readers.reader_constant import Reader as ConstantReader
from opendrift.models.larvalfish import LarvalFish

o = LarvalFish(loglevel=50)

#%%
# Seeding 20 fish eggs, which will hatch as larvae after 20 days (preliminary simplification)
time = datetime(2020, 7, 1, 12)
o.seed_elements(lon=4, lat=60, time=[time, time+timedelta(hours=24)], number=20)

#%%
# No horizontal movement, here only investigating vertical mixing and swimming
r = ConstantReader(
        {'x_sea_water_velocity': 0, 'y_sea_water_velocity': 0, 'x_wind': 0, 'y_wind': 0,
         'land_binary_mask': 0, 'ocean_vertical_diffusivity': .01})
o.add_reader(r)

o.set_config('general:use_auto_landmask', False)
o.run(duration=timedelta(days=40))

#%%
# After 20 days eggs are hatched as Larvae, and starting to grow
o.plot_property('weight')
o.plot_property('length')

#%%
# We see that larvae (after 20 days) avoid the upper meters at daytime, to avoid predators
o.plot_property('z')
