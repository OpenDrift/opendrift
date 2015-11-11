#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# HYCOM
#reader_hycom = reader_netCDF_CF_generic.Reader('http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/2010/3hrly')
#print reader_hycom

reader_globcurrent = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREUL_HS-ALT_SUM-V01.0_FULL_TIME_SERIE')

# OceanWind
reader_oceanwind = reader_netCDF_CF_generic.Reader('http://www.ncdc.noaa.gov/thredds/dodsC/oceanwinds6hr')
#print reader_oceanwind

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=100, llcrnrlat=5,
                    urcrnrlon=120, urcrnrlat=15, resolution='h')

# Add readers
print 'adding...'
o.add_reader([reader_globcurrent, reader_oceanwind, reader_basemap])


# Seed some particles
lat=10.228248; lon=106.973337
lat=10.0; lon=107.8
time = datetime(2010, 3, 23, 6, 0, 0)
o.seed_point(lon, lat, radius=1000, number=1000, time=time)

# Run model
print o
o.run(steps=200*4, time_step=900)

# Print and plot results
print o
o.plot()
