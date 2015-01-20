#!/usr/bin/env python

import datetime
import numpy as np
import matplotlib.pyplot as plt

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
#from openDriftSimulation import *
from models.od3d import OD3D

o = OD3D(proj4='+proj=stere +lat_0=90 +lon_0=70 +lat_ts=60 +units=m +a=6.371e+06 +e=0 +no_defs')


## Arome
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
#o.readers.add_reader(reader_arome, name='arome_thredds')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')
o.readers.add_reader(reader_norkyst, name='norkyst800_thredds')

lons = np.linspace(3., 3., 30)
lats = np.linspace(60., 64., 30)
depths = np.linspace(0, 0, 30)
x,y = reader_norkyst.lonlat2xy(lons, lats)
time = datetime.datetime(2015,1,15,12,0,0)
##time = datetime.datetime.now()
#
#v = reader_norkyst.get_variables(['x_sea_water_velocity',
#                                  'y_sea_water_velocity'],
#                                  time, x, y, depths, block=True)

# Arctic20
reader_arctic20 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/arctic20km/1h/aggregate_be')
o.readers.add_reader(reader_arctic20, name='arctic20_thredds')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=-3, llcrnrlat=59,
                    urcrnrlon=20, urcrnrlat=67, resolution='i')
o.readers.add_reader(reader_basemap, name='basemap_landmask')


#nork = o.readers.readers[1].get_variables(['x_sea_water_velocity'],
#                            time, x, y, depths, block=False)
#print nork

# Test Reader
#print o.readers.list_environment_variables()

#print o.proj4
#x = 10; y = 5
#print o.get_environment(['x_sea_water_velocity'], o.proj4, x, y, depths, datetime.datetime.now())
#print o.readers.get_environment(['x_sea_water_velocity'], o.proj4, x, y, depths, time)
#
o.seed_point(lon=3, lat=60, radius=10000, number=10, massOil=5, time=time)

print o.get_environment()
stop
#o.run()
#o.get_environment(['x_wind', 'y_wind', 'salinity'], 0,0,0,0)


stop


for reader in o.readers.readers:
    print reader

# Make some random points
proj4 = '+proj=latlong'
lons = np.linspace(3, 3, 10)
lats = np.linspace(60, 64, 10)
depths = np.linspace(0,0,10)
time = datetime.datetime.now()
time = datetime.datetime(2015,1,17,12,0,0)

# Norkyst
t = datetime.datetime.now()
x,y = o.readers.readers[1].lonlat2xy(lons, lats)
nork = o.readers.readers[1].get_variables(['land_binary_mask'],
                            time, x, y, depths, block=False)
print datetime.datetime.now()-t; t = datetime.datetime.now()
#
## Arctic20
#x,y = o.readers.readers[2].lonlat2xy(lons, lats)
#arc = o.readers.readers[2].get_variables(['x_sea_water_velocity'],
#                            time, x, y, depths, block=False)
#print datetime.datetime.now()-t; t = datetime.datetime.now()
#
## Basemap
#x,y = o.readers.readers[0].lonlat2xy(lons, lats)
#base = o.readers.readers[0].get_variables(['land_binary_mask'],
#                            time, x, y, depths, block=False)
#print datetime.datetime.now()-t; t = datetime.datetime.now()

print nork
print nork.shape
print x.shape
