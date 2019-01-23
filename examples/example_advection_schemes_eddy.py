#!/usr/bin/env python

# Illustrating the difference between Euler and Runge-Kutta propagation
# schemes, using an idealised (analytical) eddy current field.

# Double gyre current field from
# http://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html

import numpy as np
from datetime import datetime, timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_double_gyre
from opendrift.readers import reader_ArtificialOceanEddy
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

fake_eddy = reader_ArtificialOceanEddy.Reader(2, 62)

reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=-1.5, llcrnrlat=59,
                    urcrnrlon=7, urcrnrlat=64, resolution='i')

lon = 2.0; lat = 63.0; # Close to Station M

runs = []
leg = []
i = 0
for scheme in ['euler', 'runge-kutta', 'runge-kutta4']:
    for time_step  in [1800, 3600*3]:
        leg.append(scheme + ', T=%.1fh' % (time_step/3600.))
        print(leg[-1])
        o = OceanDrift(loglevel=50)
        o.fallback_values['land_binary_mask'] = 0
        o.set_config('drift:scheme', scheme)
        o.add_reader([fake_eddy, reader_basemap])
        o.seed_elements(lon, lat, time=datetime.now())
        o.run(duration=timedelta(days=9), time_step=time_step)
        runs.append(o)
        i = i + 1

runs[0].plot(compare=runs[1:], legend=leg)
