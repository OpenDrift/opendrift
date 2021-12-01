#!/usr/bin/env python
"""
Use Telemac (selafin) output
==================================
"""

from datetime import timedelta, datetime
from os import sep
from pyproj import Proj
from opendrift.readers import reader_telemac_selafin
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)
filename='{}Telemac_3D{}r3d_tide_open_drift.slf'.format(o.test_data_folder(),sep)
#Lambert North
proj="+proj=lcc +lat_1=49.50000000000001 +lat_0=49.50000000000001 +lon_0=0 \
     +k_0=0.999877341 +x_0=600000 +y_0=200000 +a=6378249.2 +b=6356515 \
      +units=m +no_defs"
start_time= datetime(2021,1,1,00,00)
selafin = reader_telemac_selafin.Reader(filename=filename,proj4 = proj, start_time=start_time)
o.add_reader(selafin)
o.set_config('general:coastline_action', 'previous')
# start_time = selafin.start_time generally wrong
time_step = timedelta(seconds=selafin.slf.tags["times"][1])
length=timedelta(seconds=selafin.slf.tags["times"][-1])
num_steps = len(selafin.slf.tags["times"])
# center seeds in the middle
x,y = (selafin.x.max()-selafin.x.min())/2+selafin.x.min(),\
      (selafin.y.max()-selafin.y.min())/2+selafin.y.min()
p= Proj(proj, preserve_units=False)
lon, lat = p(x,y,inverse=True)
o.seed_elements(lon=lon, lat=lat, radius=20000, number= 200, z= 0, \
   time= start_time)
o.run(time_step=time_step/10, duration=length)

o.plot(fast = True)

