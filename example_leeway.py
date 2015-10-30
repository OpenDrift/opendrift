#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.leeway import Leeway

lw = Leeway(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('/opdata_local/arome2_5/arome_metcoop_default2_5km_20150212_00.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
#reader_norkyst = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')#, name='norkyst800_file')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=-5, llcrnrlat=54,
                    urcrnrlon=27, urcrnrlat=79, resolution='h')

lw.add_reader([reader_norkyst, reader_arome, reader_basemap])

# Seeding some particles
lon = 4.5; lat = 60.0; # Outside Bergen

# Seed leeway elements at defined position and time
lw.seed_leeway(lon, lat, radius=1000, number=500,
               time=reader_arome.start_time, objectType=5)

# Running model (until end of driver data)
lw.run(steps=60*4, time_step=900, outfile='outleeway.nc')

# Print and plot results
print lw
lw.plot(buffer=.1)
