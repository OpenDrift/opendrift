#!/usr/bin/env python

# Illustrating the difference between Euler and Runge-Kutta propagation
# schemes, using an idealised (analytical) eddy current field.

from datetime import datetime

from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=50)  # Set loglevel to 0 for debug information
o.fallback_values['land_binary_mask'] = 0

double_gyre = reader_double_gyre.Reader()
print double_gyre

o.add_reader(double_gyre)


x = [.8]
y = [.5]
lon, lat = double_gyre.xy2lonlat(x, y)

o.seed_elements(lon, lat, radius=.1, number=100,
                time=double_gyre.initial_time)
o.run(steps=1000, time_step=.1)
o.plot(buffer=0)
o.animation(buffer=0)
