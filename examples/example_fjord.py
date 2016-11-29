#!/usr/bin/env python

from datetime import timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.leeway import Leeway

o = Leeway(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=5.5, llcrnrlat=61.05,
                    urcrnrlon=6.65, urcrnrlat=61.21, resolution='f',
                    projection='merc')

reader_norkyst.interpolation = 'linearND'  # Slower, but extrapolates to coast
reader_arome.interpolation = 'linearND'
o.add_reader([reader_basemap, reader_norkyst, reader_arome])

# Seed elements at defined position and time
lat = 61.117594; lon = 6.55
time = None
#time = [reader_arome.start_time,
#        reader_arome.start_time + timedelta(hours=5)]
time = reader_arome.start_time
objType = 1  # 1: Person-in-water (PIW), unknown state (mean values)
o.seed_elements(lon, lat, radius=50, number=5000, time=time, objectType=objType)

print o

# Running model (until end of driver data)
o.run(steps=66*12, time_step=300)
#stop

# Print and plot results
print o
o.plot()
o.animation()
