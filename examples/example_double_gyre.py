#!/usr/bin/env python

# Illustrating the difference between Euler and Runge-Kutta propagation
# schemes, using an idealised (analytical) eddy current field.

# Double gyre current field from
# http://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html

from datetime import datetime, timedelta

from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o.fallback_values['land_binary_mask'] = 0
o.set_config('drift:scheme', 'runge-kutta4')

double_gyre = reader_double_gyre.Reader(epsilon=.25, omega=0.628, A=0.1)
print double_gyre

o.add_reader(double_gyre)

x = [.9]
y = [.5]
lon, lat = double_gyre.xy2lonlat(x, y)

o.seed_elements(lon, lat, radius=.1, number=1000,
                time=double_gyre.initial_time)

o.run(duration=timedelta(seconds=10), time_step=0.1)
o.animation(buffer=0)
