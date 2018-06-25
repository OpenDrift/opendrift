#!/usr/bin/env python

from datetime import timedelta

import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=4.0, llcrnrlat=59.9,
                    urcrnrlon=5.5, urcrnrlat=61.2,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_norkyst])

# Seeding some particles
lons = np.linspace(4.4, 4.6, 10)
lats = np.linspace(60.0, 60.1, 10)
lons, lats = np.meshgrid(lons, lats)
lons = lons.ravel()
lats = lats.ravel()

# Seed oil elements on a grid at regular time interval
start_time = reader_norkyst.start_time
time_step = timedelta(hours=6)
num_steps = 10
for i in range(num_steps+1):
    o.seed_elements(lons, lats, radius=0, number=100,
                    time=start_time + i*time_step)

# Running model (until end of driver data)
o.run(steps=66*4, time_step=900)

# Print and plot results
print(o)
o.animation()
