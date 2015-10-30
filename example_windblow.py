#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic

from models.windblow import WindBlow

o = WindBlow(loglevel=0)  # Set loglevel to 0 for debug information

reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('/opdata_local/arome2_5/arome_metcoop_default2_5km_20150727_18.nc')

#print reader_arome
#stop

reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=4.5, llcrnrlat=59,
                    urcrnrlon=12, urcrnrlat=65, resolution='h')

o.add_reader([reader_arome, reader_basemap])

# Seeding some particles
lon = 4.5; lat = 63.0; 
lon = 9.6; lat = 64.3; 
o.seed_point(lon, lat, radius=50000, number=500, time=reader_arome.start_time)

# Running model (until end of driver data)
o.run(steps=125*4, time_step=900)

# Print and plot results
print o
o.plot(buffer=.1)
