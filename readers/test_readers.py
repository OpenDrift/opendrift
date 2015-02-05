#!/usr/bin/env python

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

import reader_basemap_landmask
import reader_netCDF_CF_generic
from readers import Readers

r = Readers()

# Arome
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')  #, name='arome_thredds')
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

r.add_reader([reader_norkyst, reader_arctic20], ['x_sea_water_velocity', 'y_sea_water_velocity'])
r.add_reader([reader_arome])
r.add_reader([reader_arctic20])
r.add_reader([reader_norkyst, reader_basemap])
#r.add_reader([reader_arome], ['low_type_cloud_area_fraction', 'air_temperature'])
print r

#print r.readers['http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc'].variables
