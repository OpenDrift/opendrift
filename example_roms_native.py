#!/usr/bin/env python

import os
from datetime import datetime, timedelta
import numpy as np

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from readers import reader_ROMS_native
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Nordic 4km
file = '/disk2/data/roms_curvilinear/tilKnutFrode/ocean_his.nc'
if not os.path.exists(file):
    raise ValueError('Native ROMS file not available: ' + file)
nordic_native = reader_ROMS_native.Reader(file)

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=10.3, llcrnrlat=59.0,
                    urcrnrlon=11.0, urcrnrlat=59.7,
                    resolution='f', projection='merc')

o.add_reader([reader_basemap, nordic_native])

# Seeding some particles
time = nordic_native.start_time
#lons, lats = nordic_native.domain_grid()
#lon = lons.ravel()
#lat = lats.ravel()
lon = 10.52; lat = 59.46; # Oslofjord

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=2000, number=5000, time=time)

print o

# Adjusting some configuration
o.config['drift']['wind_drift_factor'] = .02
o.config['processes']['diffusion'] = True
o.config['processes']['dispersion'] = True
o.config['processes']['evaporation'] = True
o.config['processes']['emulsification'] = True
o.config['drift']['current_uncertainty'] = .1
o.config['drift']['wind_uncertainty'] = 2

# Running model
o.run(steps=24*4*3, time_step=900)

# Print and plot results
print o
o.animation()
o.plot()
