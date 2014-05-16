#!/usr/bin/env python

import numpy as np
import reader_norkystOpenDap

r = reader_norkystOpenDap.Reader()

lon = np.array([3, 3.2, 3])
lat = np.array([60, 60.1, 60.0])
x, y = r.lonlat2xy(lon, lat)
print x, y
res = r.get_parameters('sea_water_temperature', time=None,
                        x=x, y=y, depth=None)
print res

print r
