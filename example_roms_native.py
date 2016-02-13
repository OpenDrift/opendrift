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
nordic_native = reader_ROMS_native.Reader('test_data/2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=10.0, llcrnrlat=67.5,
                    urcrnrlon=16.0, urcrnrlat=69.0,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, nordic_native])

# Seeding some particles
time = nordic_native.start_time
lon = 12.0; lat = 68.3;

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
o.run(steps=24*2, time_step=3600)

# Print and plot results
print o
o.animation()
o.plot()
