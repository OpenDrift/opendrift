#!/usr/bin/env python

from datetime import datetime
import numpy as np
import reader_norkystOpenDap
import reader_ROMSArctic20openDap
import reader_netCDF_CF_generic

n800 = reader_norkystOpenDap.Reader()
a20 = reader_ROMSArctic20openDap.Reader()
ng = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')
#ng = reader_netCDF_CF_generic.Reader('/opdata/arome_norway25/AROME_Norway25_00.nc')

r = ng
r = n800
r = a20
print r
#
lon = np.array([1, 3.2, 9.50])
lat = np.array([60, 59.1, 65.042])
t = datetime(2014, 5, 22, 2, 0, 0)
x, y = r.lonlat2xy(lon, lat)
res = r.get_parameters('sea_water_salinity', time=t,
                        x=x, y=y, depth=None)#
print r
print res
#indx, time = r.index_of_closest_time(t)
