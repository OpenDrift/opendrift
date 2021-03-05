#!/usr/bin/env python
"""
Double gyre
=============

Drift in an idealised (analytical) double gyre eddy current field from
https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html
"""

from datetime import datetime, timedelta

from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o.set_config('environment:fallback:land_binary_mask', 0)
o.set_config('drift:advection_scheme', 'runge-kutta4')

double_gyre = reader_double_gyre.Reader(epsilon=.25, omega=0.628, A=0.1)
print (double_gyre)

o.add_reader(double_gyre)

x = [.9]
y = [.5]
lon, lat = double_gyre.xy2lonlat(x, y)

o.seed_elements(lon, lat, radius=.1, number=5000,
                time=double_gyre.initial_time)

o.run(duration=timedelta(seconds=10), time_step=0.1)
o.animation(buffer=0, hide_landmask=True)

#%%
# .. image:: /gallery/animations/example_double_gyre_0.gif

o.plot(buffer=0, hide_landmask=True)
