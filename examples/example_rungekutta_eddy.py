#!/usr/bin/env python

# Illustrating the difference between Euler and Runge-Kutta propagation
# schemes, using an idealised (analytical) eddy current field.

from datetime import datetime

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_ArtificialOceanEddy
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

fake_eddy = reader_ArtificialOceanEddy.Reader(2, 62)

reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=-1.5, llcrnrlat=59,
                    urcrnrlon=7, urcrnrlat=64, resolution='h')

o.add_reader([fake_eddy, reader_basemap])
lon = 2.0; lat = 63.0; # Close to Station M

# First run, with Euler scheme:
o.set_config('drift:scheme', 'euler')
o.seed_elements(lon, lat, radius=0, number=1, time=datetime(2015,1,1))
o.run(steps=300, time_step=3600)

# Second run, with Runge-Kutta scheme:
o2 = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o2.add_reader([fake_eddy, reader_basemap])
o2.set_config('drift:scheme', 'runge-kutta')
o2.seed_elements(lon, lat, radius=0, number=1, time=datetime(2015,1,1))
o2.run(steps=300, time_step=3600)

o.plot(compare=o2, legend=['Euler scheme', 'Runge-Kutta scheme'])
#o2.plot()
# Animate and compare the two runs
o.animation(compare=o2, legend=['Euler scheme', 'Runge-Kutta scheme'])
