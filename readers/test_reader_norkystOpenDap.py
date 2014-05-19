#!/usr/bin/env python

from datetime import datetime
import numpy as np
import reader_norkystOpenDap
import reader_ROMSArctic20openDap

n800 = reader_norkystOpenDap.Reader()
a20 = reader_ROMSArctic20openDap.Reader()

r = a20

lon = np.array([3, 3.2, 3])
lat = np.array([60, 60.1, 89.0])
t = datetime(2014, 5, 19, 10, 20, 0)
x, y = r.lonlat2xy(lon, lat)
res = r.get_parameters('x_sea_water_velocity', time=t,
                        x=x, y=y, depth=None)
#print r
print res
indx, time = r.index_of_closest_time(t)
