#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic

from models.oceandrift import OceanDrift

ncfile = 'backandforth.nc'

back = False
#back = True

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')

reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=2, llcrnrlat=58,
                urcrnrlon=8, urcrnrlat=63, resolution='h')

o.add_reader([reader_norkyst, reader_arome, reader_basemap])


if back:
    o.io_import_file(ncfile)
    o.start_time = o.time  # use last time as start now

    # Running model backwards from end of forward simulation
    o.run(steps=70*4, time_step=-900, outfile='back.nc')
    figfile = 'back.png'

else:
    # Seeding some particles
    lon = 4.4; lat = 60.0; 
#    lon = 3.4; lat = 60.0; 
    time = reader_norkyst.start_time
    o.seed_point(lon, lat, radius=10000, number=100, time=time)

    # Running model backwards from end of forward simulation
    o.run(steps=70*4, time_step=900, outfile=ncfile)
    figfile = 'forward.png'

# Print and plot results
print o

o.plot(buffer=.2, filename=figfile)
