#!/usr/bin/env python

from datetime import timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.leeway import Leeway

o = Leeway(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

# Making customised, full resolution landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=5.5, llcrnrlat=61.05,
                    urcrnrlon=6.65, urcrnrlat=61.21, resolution='f',
                    projection='merc')

o.add_reader([reader_basemap, reader_norkyst, reader_arome])

# Seed elements
lat = 61.117594; lon = 6.55
#time = [reader_arome.start_time,
#        reader_arome.start_time + timedelta(hours=5)]
time = reader_arome.start_time
objType = 1  # 1: Person-in-water (PIW), unknown state (mean values)
o.seed_elements(lon, lat, radius=50, number=5000, time=time, objectType=objType)

# Running model for 66 hours
o.run(steps=66*12, time_step=300)

# Print and plot results
print(o)
o.plot()
o.animation()
