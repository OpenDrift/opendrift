#!/usr/bin/env python

import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Nordic 4km
nordic_native = reader_ROMS_native.Reader(o.test_data_folder() +
    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=11.0, llcrnrlat=67.5,
                    urcrnrlon=16.0, urcrnrlat=69.0,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, nordic_native])

# Seeding some particles
time = nordic_native.start_time
lon = 12.0; lat = 68.3;

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=0, number=10, z=np.linspace(0, -150, 10), time=time)

print(o)

# Running model
o.run(time_step=3600)

# Print and plot results
print(o)
o.plot(linecolor='z')
#o.animation()
