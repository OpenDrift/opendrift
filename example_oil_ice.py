#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')

# Arctic20
reader_arctic20 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/arctic20km/1h/aggregate_be')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=2, llcrnrlat=78,
                    urcrnrlon=22, urcrnrlat=82.3,
                    resolution='i', projection='merc')

o.add_reader([reader_basemap, reader_arctic20])
o.fallback_values['y_wind'] = 35  # Northerly wind, to get stranding in ice

# Seeding some particles
lon = 8; lat = 79.5;  # West of Svalbard

#time = datetime(2015, 9, 22, 6, 0, 0)
#time = [reader_arctic20.start_time,
#        reader_arctic20.start_time + timedelta(hours=30)]
time = reader_arctic20.start_time

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=5000, number=3000, time=time)

print o

# Adjusting some configuration
o.config['drift']['wind_drift_factor'] = .02
o.config['processes']['diffusion'] = True
o.config['processes']['dispersion'] = False
o.config['processes']['evaporation'] = False
o.config['processes']['emulsification'] = False
o.config['drift']['current_uncertainty'] = .5
o.config['drift']['wind_uncertainty'] = 5

# Running model (until end of driver data)
o.run(steps=120, time_step=3600)

# Print and plot results
print o
o.animation()
o.plot(background='sea_ice_area_fraction')
