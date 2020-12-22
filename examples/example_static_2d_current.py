#!/usr/bin/env python
"""
Static 2D current field
=======================
"""

from datetime import datetime, timedelta
import numpy as np

from opendrift.readers.reader_constant_2d import Reader
from opendrift.models.oceandrift import OceanDrift

#%%
# Constructing a static, rotating ocean current field, 
lon, lat = np.meshgrid(np.linspace(2,6,30), np.linspace(59,62,30))
lon0 = 4
lat0 = 60.5
u = -(lat-lat0)/np.sqrt((lon-lon0)**2 + (lat-lat0)**2)
v = (lon-lon0)/np.sqrt((lon-lon0)**2 + (lat-lat0)**2)
lon = np.linspace(0,5,30)
lat = np.linspace(59,62,30)

r = Reader(x=lon, y=lat, proj4='+proj=latlong',
           array_dict = {'x_sea_water_velocity': u, 'y_sea_water_velocity': v})

o = OceanDrift(loglevel=20)
o.set_config('environment:fallback:land_binary_mask', 0)
o.add_reader(r)
o.seed_elements(lon=3, lat=60.5, number=1000, radius=30000, time=datetime.now())
o.run(duration=timedelta(hours=72))
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_static_2d_current_0.gif

