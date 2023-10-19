#!/usr/bin/env python
"""
Advection schemes
=================
Illustrating the difference between Euler and Runge-Kutta propagation
schemes, using an idealised (analytical) eddy current field.
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_ArtificialOceanEddy
from opendrift.models.oceandrift import OceanDrift

fake_eddy = reader_ArtificialOceanEddy.Reader(2, 62)

runs = []
leg = []
for scheme in ['euler', 'runge-kutta', 'runge-kutta4']:
    for time_step  in [1800, 3600*3]:
        leg.append(scheme + ', T=%.1fh' % (time_step/3600.))
        print(leg[-1])
        o = OceanDrift(loglevel=50)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('drift:advection_scheme', scheme)
        o.set_config('drift:vertical_mixing', False)
        o.add_reader(fake_eddy)
        o.seed_elements(lon=2.0, lat=63.0, time=datetime.utcnow())
        o.run(duration=timedelta(days=9), time_step=time_step)
        runs.append(o)

runs[0].plot(compare=runs[1:], legend=leg, fast=True, buffer=.3)
