#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from readers import reader_ROMS_native
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Nordic 4km
#nordic_native = reader_ROMS_native.Reader('/disk2/data/roms/ocean_his_0170.nc')
nordic_native = reader_ROMS_native.Reader('/opdata/roms/Nordic-4km_SLEVELS_avg_00.nc')
#nordic_cf = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/fou-hi/nordic4km-1h/Nordic-4km_SURF_1h_avg_00.nc')

# Landmask (Basemap)
#reader_basemap = reader_basemap_landmask.Reader(
#                    llcrnrlon=13.4, llcrnrlat=67.8,
#                    urcrnrlon=14.9, urcrnrlat=68.2,
#                    resolution='h', projection='merc')
#reader_basemap = reader_basemap_landmask.Reader(
#                    llcrnrlon=3, llcrnrlat=59,
#                    urcrnrlon=9, urcrnrlat=64.0,
#                    resolution='h', projection='merc')
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=5.9, llcrnrlat=63.9,
                    urcrnrlon=16.2, urcrnrlat=68.6,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, nordic_native])

# Seeding some particles
#lon = 14.1; lat = 68.04; # Vestfjorden
#lon = 4.8; lat = 60.0; # Bergen
#lon = 10.3; lat = 64.95; # Nordland
time = nordic_native.start_time
lons = np.linspace(9, 11.0, 50)
lats = np.linspace(64, 66.5, 75)
#lons = np.linspace(3, 5.0, 50)
#lats = np.linspace(59, 61, 75)
lons, lats = np.meshgrid(lons, lats)
lon = lons.ravel()
lat = lats.ravel()

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=0, number=5000, time=time)

print o

# Adjusting some configuration
o.config['drift']['wind_drift_factor'] = .02
o.config['processes']['diffusion'] = False
o.config['processes']['dispersion'] = True
o.config['processes']['evaporation'] = True
o.config['processes']['emulsification'] = True
o.config['drift']['current_uncertainty'] = .1
o.config['drift']['wind_uncertainty'] = 2

# Running model (until end of driver data)
#o.run(steps=20, time_step=1800, outfile='native.nc')
o.run(steps=24*4, time_step=900)

# Print and plot results
print o
o.animation()
#o.plot()
