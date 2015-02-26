#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from readers import reader_netcdf_new
from models.oceandrift import OceanDrift
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# HYCOM
reader_hycom = reader_netcdf_new.Reader('http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/2010/3hrly')
print reader_hycom

reader_oceanwind = reader_netcdf_new.Reader('http://www.ncdc.noaa.gov/thredds/dodsC/oceanwinds6hr')
print reader_oceanwind

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=-95, llcrnrlat=20,
                    urcrnrlon=-70, urcrnrlat=40, resolution='i')

# Add readers
print 'adding...'
o.add_reader([reader_hycom, reader_oceanwind, reader_basemap])

# Seed some particles
lon = -88.387161; lat = 28.736669  # Macondo location
time = datetime(2010, 4, 21, 6, 0, 0)  # 4 hours after explosion
o.seed_point(lon, lat, radius=10000, number=500, massOil=1, time=time)

# Run model
print o
o.run(steps=20, time_step=3600*3)  # Using three hour timestep of hycom

# Print and plot results
print o
#o.plot(background='sea_water_potential_temperature')
o.plot()
