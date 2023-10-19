#!/usr/bin/env python
"""
Sediment drift with resuspension
================================
"""

from datetime import timedelta, datetime
from opendrift.readers import reader_oscillating
from opendrift.models.sedimentdrift import SedimentDrift

#%%
# Constructing an artificial current field where x- and y-components are oscilating with different amplitude and period
reader_oscx = reader_oscillating.Reader('x_sea_water_velocity',
    amplitude=0.6, zero_time=datetime.utcnow())
reader_oscy = reader_oscillating.Reader('y_sea_water_velocity',
    amplitude=.3, period_seconds=3600*5, zero_time=datetime.utcnow())

o = SedimentDrift(loglevel=50)  # 0 for debug output

if True:
    o.add_reader([reader_oscx, reader_oscy])
    o.set_config('environment:fallback:y_wind', -6)
    o.set_config('environment:fallback:x_wind', -3)
    o.set_config('environment:fallback:sea_floor_depth_below_sea_level', 30)  # 100m depth
else:  # Using live data from Thredds instead of oscillating currents
    o.add_readers_from_list([
        'https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be'])

#%%
# Set threshold for bottom resuspension
o.set_config('vertical_mixing:resuspension_threshold', .5)

# Adding some horizontal and vertical diffusion
o.set_config('drift:current_uncertainty', 0.1)
o.set_config('drift:wind_uncertainty', 1)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Large1994')
#o.set_config('vertical_mixing:diffusivitymodel', 'environment')

#%%
# Seeding sediments
o.seed_elements(lon=4.65, lat=60, number=10000, 
                time=[datetime.utcnow(), datetime.utcnow()+timedelta(hours=6)],
                terminal_velocity=-.01)  # 1 cm/s settling speed

o.run(time_step=1800, time_step_output=1800, duration=timedelta(hours=72))

#%%
# Plotting the depth vs time
o.plot_property('z')

#%%
# Animate sediment particles
o.animation(color='moving', colorbar=False, legend=['Sedimented', 'Moving'], fast=True, buffer=.01)
#o.animation_profile()

#%%
# .. image:: /gallery/animations/example_sediments_resuspension_0.gif
