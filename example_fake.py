#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_ArtificialOceanEddy
from models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0, time_step=300)  # Set loglevel to 0 for debug information

fake_eddy = reader_ArtificialOceanEddy.Reader(2, 62)
#fake_eddy.plot()
o.add_reader([fake_eddy])


reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=0, llcrnrlat=57,
                    urcrnrlon=10, urcrnrlat=67, resolution='h')
o.add_reader([reader_basemap])
print reader_basemap

print o

# Seeding some particles
lon = 2.0; lat = 65.0; # Close to Station M
lon = 3.0; lat = 60.0; # Off Bergen, for stranding test
o.seed_point(lon, lat, radius=10000, number=100, time=None)

# Running model (until end of driver data)
o.run(steps=6000)

# Print and plot results
print o
o.plot(buffer=2)
