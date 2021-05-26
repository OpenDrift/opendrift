# tuto
from datetime import timedelta
import numpy as np
from os import sep
from pyproj import Proj
from opendrift.readers import reader_telemac_selafin
from opendrift.readers import reader_global_landmask
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)
filename='{}Telemac_3D{}r3d_tide_open_drift.slf'.format(o.test_data_folder(),sep)
#Lambert North
proj="+proj=lcc +lat_1=49.50000000000001 +lat_0=49.50000000000001 +lon_0=0 \
     +k_0=0.999877341 +x_0=600000 +y_0=200000 +a=6378249.2 +b=6356515 \
      +units=m +no_defs"
selafin = reader_telemac_selafin.Reader(filename=filename,proj4 = proj)
o.add_reader(selafin)
o.set_config('general:coastline_action', 'previous')
start_time = selafin.start_time
time_step = timedelta(seconds=selafin.slf.tags["times"][1])
length=timedelta(seconds=selafin.slf.tags["times"][-1])
num_steps = len(selafin.slf.tags["times"])

# center seeds in the middle
lon, lat = selafin.center()

o.seed_elements(lon=lon, lat=lat, radius=20000, number= 200, z= 0, \
   time= start_time)
o.run(time_step=time_step/10, duration=length)
