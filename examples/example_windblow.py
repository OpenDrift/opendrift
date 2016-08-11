#!/usr/bin/env python

from datetime import datetime, timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic

from opendrift.models.windblow import WindBlow

o = WindBlow(loglevel=0)  # Set loglevel to 0 for debug information

#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=2.5, llcrnrlat=59,
                    urcrnrlon=8, urcrnrlat=63, resolution='h')

o.add_reader([reader_arome, reader_basemap])

# Seeding some particles
lon = 6.5; lat = 60.3; 
o.seed_elements(lon, lat, radius=5000, number=1000,
                time=reader_arome.start_time)

# Running model (until end of driver data)
o.run(steps=66*4, time_step=900)

# Print and plot results
print o
o.plot(buffer=.5)
