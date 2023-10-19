#!/usr/bin/env python
"""
Comparing Leeway and ShipDrift model
====================================
"""

#%%
# The Leeway model contain many object categories from persons-in-water
# to various types of boats and ships.
# Here we compare the Leeway model with a fishing vessel to
# the ShipDrift model using the same ship dimensions.

from datetime import datetime, timedelta
from opendrift.models.leeway import Leeway
from opendrift.models.shipdrift import ShipDrift

#%%
# We use a simple case with constant wind northwards, and no current.
# Wave height and period is calculated automatically from wind,
# and wave direction is the same as wind direction
environment = {
    'land_binary_mask': 0,
    'x_sea_water_velocity': 0,
    'y_sea_water_velocity': 0,
    'x_wind': 0,
    'y_wind': 18}

#%%
# For the Leeway model we use "FISHING-VESSEL-3" with a length of 62m
# For the ShipDrift model we define a ship with the same dimensions,
# as well as a "large ship" scaled up by a factor of 8
ship = {'length': 62, 'beam': 8, 'height': 10, 'draft': 5}
large_ship = {k: v*8 for k, v in ship.items()}

cases = {
    'Leeway, FISHING-VESSEL-3': {
        'model': Leeway, 'kwargs': {'object_type': 52}},
    'Shipdrift, large ship': {
        'model': ShipDrift, 'kwargs': large_ship},
    'Shipdrift, small ship': {
        'model': ShipDrift, 'kwargs': ship}
    }

lon=3.5
lat=60
time = datetime.now()
duration = timedelta(hours=24)
simulations = []
for cname, case in cases.items():
    o = case['model'](loglevel=50)
    for var, value in environment.items():
        o.set_config('environment:constant:' + var, value)
    o.set_config('general:use_auto_landmask', False)
    o.seed_elements(lon=lon, lat=lat, time=time, number=1000,
                    **case['kwargs'])
    o.run(duration=duration)
    simulations.append(o)

#%%
# We see that the overall drift speed of the ShipDrift model
# is comparable to the Leeway model.
# However, the Leeway model yields a much larger directional
# spread for ships orienters left and right respectively.
# For the scaled up ship, the drift speed is slightly larger,
# and the directional spread is slightly larger, but still
# much smaller than with the Leeway model.
simulations[0].plot(compare=simulations[1:], legend=list(cases), fast=False)
simulations[0].animation(compare=simulations[1:], legend=list(cases), fast=False)

#%%
# .. image:: /gallery/animations/example_shipdrift_leeway_0.gif
