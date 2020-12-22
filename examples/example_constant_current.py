#!/usr/bin/env python
"""
Constant current
================
"""

from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Adding no input models, but instead constant northwards current of 1 m/s
o.set_config('environment:fallback:x_sea_water_velocity', 0)
o.set_config('environment:fallback:y_sea_water_velocity', 1)
o.set_config('environment:fallback:land_binary_mask', 0)

#%%
# Seed elements at defined position and time
o.seed_elements(lon=4.0, lat=60.0, radius=5000, number=100,
                time=datetime(2015, 9, 22, 6, 0, 0))

#%%
# Running model for 50 hours
o.run(duration=timedelta(hours=50))

#%%
# Print and plot results
print(o)
o.plot(fast=True, buffer=1)
