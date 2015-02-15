#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil
from models.windblow import WindBlow

#o = WindBlow()
#now = datetime.now()
#o.seed_point(4.9, 60.0, radius=1000, number=10, time=now, length=1, CodLarvaeProperty1=1)
#print o.elements
##o.seed_point(4.9, 60.0, radius=1000, number=10, time=now)
##print o.elements
#stop

#o = OpenDriftSimulation()
o = OpenOil()

#from test_concrete import *
#c = Concrete()
#c.foo()
#c.bar()
#stop

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')  #, name='arome_thredds')
#reader_arome = reader_netCDF_CF_generic.Reader('/opdata_local/arome2_5/arome_metcoop_default2_5km_20150212_00.nc')
#o.readers.add_reader(reader_arome)

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
#reader_norkyst = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')#, name='norkyst800_file')

# Arctic20
reader_arctic20 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/arctic20km/1h/aggregate_be', name='arctic20_thredds')
#o.add_reader(reader_arctic20)


# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=-5, llcrnrlat=54,
                    urcrnrlon=20, urcrnrlat=69, resolution='i')
#o.add_reader(reader_basemap)

reader_basemap.plot()
reader_norkyst.plot()
reader_arctic20.plot()

#o.add_reader(reader_arome)
#o.add_reader([reader_norkyst, reader_arctic20], ['x_sea_water_velocity', 'y_sea_water_velocity'])
o.add_reader([reader_norkyst])

print o

# Seeding some particles
#lon = 15; lat = 72.0; # Close to Norkyst boundary
#lon = 21; lat = 73.5; # Close to Norkyst boundary
lon = 4.9; lat = 60.0; # Outside Bergen
o.seed_point(lon, lat, radius=10000, number=5, massOil=5, time=None)

# Running model (until end of driver data)
o.use_block = True
#o.time_step = timedelta(seconds=900)
o.run(steps=100)

# Print and plot results
print o
#o.plot(background='sea_water_potential_temperature')
o.plot()
