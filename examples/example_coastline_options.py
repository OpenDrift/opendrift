#!/usr/bin/env python
"""
Coastline interaction
======================

Example to illustrate stranding options using an artificial
east-west oscillating current field
Knut-Frode Dagestad, Feb 2017
"""

from datetime import datetime
from opendrift.readers import reader_oscillating
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=50)  # Set loglevel to 0 for debug information

reader_osc = reader_oscillating.Reader('x_sea_water_velocity', amplitude=1,
                                       zero_time=datetime.utcnow())
o.add_reader([reader_osc])  # Oscillating east-west current component

o.fallback_values['y_sea_water_velocity'] = .2  # Adding northwards drift

#%%
# Try different options: 'previous', 'stranding', 'none'
o.set_config('general:coastline_action', 'previous')

o.seed_elements(lon=12.2, lat=67.7, radius=5000, number=25, time=reader_osc.zero_time)

o.run(steps=36*4, time_step=900, time_step_output=1800)

#%%
# Print and plot results
print(o)
o.animation()
o.plot()

#%%
# .. image:: /gallery/animations/example_coastline_options_0.gif
