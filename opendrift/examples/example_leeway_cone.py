#!/usr/bin/env python

from datetime import datetime, timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.leeway import Leeway

lw = Leeway(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
reader_arome = reader_netCDF_CF_generic.Reader(
    '../test_data/16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst = reader_netCDF_CF_generic.Reader(
    '../test_data/16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=2.5, llcrnrlat=59.3,
                    urcrnrlon=5.8, urcrnrlat=62.5, resolution='h')

lw.add_reader([reader_norkyst, reader_arome, reader_basemap])

# Seed elements along cone, e.g. ship track with
# increasing uncertainty in position
lon = [3.6, 5.1]; lat = [61., 59.6];
time = [reader_arome.start_time, reader_arome.start_time + timedelta(hours=30)]
#time = reader_arome.start_time

objType = 26  # 26 = Life-raft, no ballast
lw.seed_elements(lon, lat, radius=[1000, 10000], number=5000,
                 time=time, objectType=objType)

# Running model (until end of driver data)
lw.run(steps=66*4, time_step=900, outfile='outleeway.nc')

# Print and plot results
print lw
lw.plot()
lw.animation()
