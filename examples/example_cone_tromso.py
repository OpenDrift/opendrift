#!/usr/bin/env python

from datetime import timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('/disk1/data/opendrift_testdata/arome.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
#reader_norkyst = reader_netCDF_CF_generic.Reader('/disk1/data/opendrift_testdata/norkyst800.nc')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=15.4, llcrnrlat=68.6,
                    urcrnrlon=20.0, urcrnrlat=70.5,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_norkyst, reader_arome])

# Seeding some particles

latstart = 68.988911 
lonstart = 16.040701
#latstart = 69.477754
#lonstart = 16.441702
latend = 69.991446
lonend = 17.760061
lon = [lonstart, lonend]; lat = [latstart, latend]; # Outside Tromso

# Seed elements along cone, e.g. ship track with
# increasing uncertainty in position
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]
#time = reader_arome.start_time

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=[100, 500], number=10000,
                time=time, cone=True)

print o

# Adjusting some configuration
o.config['processes']['diffusion'] = True
o.config['processes']['dispersion'] = True
o.config['processes']['evaporation'] = False
o.config['processes']['emulsification'] = True

# Running model (until end of driver data)
o.run(steps=66*2, time_step=1800)

# Print and plot results
print o
o.plot()
o.animation()
#o.animation(filename='oilspill_tromsoe.mp4')
