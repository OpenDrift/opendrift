#!/usr/bin/env python
"""
Sediment drift from Glomma river outlet
=======================================
"""

from datetime import timedelta, datetime
from opendrift.models.sedimentdrift import SedimentDrift

o = SedimentDrift(loglevel=50)  # 0 for debug output

#%%
# Source of sediments at Glomma river outlet
lat=59.169194
lon=10.962920
o.seed_elements(lon=lon, lat=lat, number=10000, 
                time=[datetime.utcnow(), datetime.utcnow()+timedelta(hours=48)],
                terminal_velocity=-.001)  # 1 mm/s settling speed

if False:  # Using constant south-westwards current and wind
    o.fallback_values['x_sea_water_velocity'] = -.05
    o.fallback_values['y_sea_water_velocity'] = -.1
    o.fallback_values['y_wind'] = -6
    o.fallback_values['x_wind'] = -2
    o.fallback_values['sea_floor_depth_below_sea_level'] = 100  # 100m depth
else:  # Using live data from Thredds
    o.add_readers_from_list([
        'https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be'])

#%%
# Adding some diffusion
o.set_config('drift:current_uncertainty', .2)
o.set_config('drift:wind_uncertainty', 2)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Large1994')
#o.set_config('vertical_mixing:diffusivitymodel', 'environment')

o.run(time_step=600, time_step_output=1800, duration=timedelta(hours=48))

#%%
# Plotting the depth vs time
o.plot_property('z')

#%%
# Animate sediment particles colored by their depth
o.animation(color='z', fast=False, buffer=.01)
o.animation(color='moving', fast=False, buffer=.01)

#%%
# .. image:: /gallery/animations/example_sediments_0.gif
