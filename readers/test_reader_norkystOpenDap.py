#!/usr/bin/env python

import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import reader_norkystOpenDap
import reader_ROMSArctic20openDap
import reader_netCDF_CF_generic
import reader_basemap_landmask

n800_od = reader_norkystOpenDap.Reader()
a20 = reader_ROMSArctic20openDap.Reader()
ng = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')
#ng = reader_netCDF_CF_generic.Reader('/opdata/arome_norway25/AROME_Norway25_00.nc')
r = reader_basemap_landmask.Reader(llcrnrlon=-3, llcrnrlat=59,
                                   urcrnrlon=10, urcrnrlat=67, resolution='i')
r = ng
#r = n800_od
#r = a20
print r

lon = np.array([3, 3.2, 5.50])
lat = np.array([60, 59.1, 65.042])
#size = 10000
#lon = np.random.uniform(-3, 10, size)
#lat = np.random.uniform(59, 67, size)

t = datetime(2014, 11, 19, 2, 0, 0)
t = datetime.now()
x, y = r.lonlat2xy(lon, lat)

res = a20.get_parameters('sea_water_salinity', time=t, x=x, y=y, depth=None)#
res2 = n800_od.get_parameters('sea_water_salinity', time=t, x=x, y=y, depth=None)#
print r
print res
print res2

r.map.drawcoastlines()
r.map.drawcountries(linewidth=0.25)
r.map.fillcontinents(color='coral',lake_color='aqua')
r.map.drawmeridians(np.arange(0,360,3))
r.map.drawparallels(np.arange(-90,90,3))
#r.map.plot(x, y, '*')
r.map.plot(lon[res], lat[res], 'r.')
r.map.plot(lon[~res], lat[~res], 'b.')
plt.show()
