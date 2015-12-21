#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from readers import reader_ROMS_native
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Nordic 4km
#reader_nordic4 = reader_netCDF_CF_generic.Reader('test_data/norkyst800_subset_16Nov2015.nc')
reader_nordic4 = reader_ROMS_native.Reader('/disk2/data/roms/ocean_his_0170.nc')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=13.4, llcrnrlat=67.8,
                    urcrnrlon=14.9, urcrnrlat=68.2,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_nordic4])

# Seeding some particles
lon = 14.1; lat = 68.04; # Vestfjorden
time = reader_nordic4.start_time

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, number=4, time=time)

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
o.run(steps=24*2, time_step=1800, outfile='roms.nc')

# Print and plot results
print o
o.plot()
#o.animation()
