#!/usr/bin/env python

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.od3d import OD3D
from models.windblow import WindBlow

#o = OD3D(proj4='+proj=stere +lat_0=90 +lon_0=70 +lat_ts=60 +units=m +a=6.371e+06 +e=0 +no_defs')
#o = OD3D()
o = WindBlow()

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')  #, name='arome_thredds')
reader_arome = reader_netCDF_CF_generic.Reader('/opdata_local/arome2_5/arome_metcoop_default2_5km_20150205_00.nc')
#o.readers.add_reader(reader_arome)

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')#, name='norkyst800_file')

# Arctic20
reader_arctic20 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/arctic20km/1h/aggregate_be', name='arctic20_thredds')
#o.add_reader(reader_arctic20)

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=-5, llcrnrlat=54,
                    urcrnrlon=20, urcrnrlat=69, resolution='i')
#o.add_reader(reader_basemap)

#o.add_reader([reader_norkyst, reader_arctic20], ['x_sea_water_velocity', 'y_sea_water_velocity'])
#o.add_reader([reader_norkyst, reader_arctic20, reader_arome, reader_basemap])
o.add_reader(reader_arome)

# Seeding some particles
#time = datetime(2015,1,20,1,0,0) # Arctic20
#time = datetime(2015,1,15,0,0,0) # Norkyst800
#o.seed_point(lon=4.8, lat=60.0, radius=10000, number=30, massOil=5, time=None)
o.seed_point(lon=8, lat=60.0, radius=20000, number=30, time=None)
print o

# Running model (until end of driver data)
o.run()

# Print and plot results
print o
#o.plot('x_sea_water_velocity')
o.plot('sea_water_potential_temperature')
