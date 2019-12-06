#!/usr/bin/env python
"""
CMEMS
=============

This example runs a Leeway simulation, with current data downloaded
from CMEMS, and wind data from an NCEP Thredds server.
To run this example, you need a CMEMS account created at
http://marine.copernicus.eu, and you need to install the motuclient
available from https://github.com/clstoulouse/motu-client-python :
python -m pip install motuclient
"""

import os
from datetime import datetime, timedelta
import numpy as np

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_cmems
from opendrift.readers import reader_netCDF_CF_generic

try:
    f = open('cmems_account.txt')
except:
    raise ValueError('Please store your CMEMS username and password '
                     'in local file "cmems_account.txt"')

username = f.readline().strip()
password = f.readline().strip()
f.close()

# Seed information
lon = 4.8; lat = 60  # Bergen, Norway
lon = -89; lat = 29.8  # New Orleans
lon = 107; lat = 10  # Ho Chi Minh

time = datetime.now()
duration = timedelta(days=3)
bufferlat = duration.total_seconds()/111000
bufferlon = bufferlat*np.cos(lat*np.pi/180)

# Fetching current data from CMEMS
cmems_file = 'opendrift_cmems_download.nc'
if os.path.exists(cmems_file):
    # Reuising downloaded file, if existing. Delete it to force update.
    cmems = reader_netCDF_CF_generic.Reader(cmems_file)
else:
    cmems = reader_cmems.Reader(username=username, password=password,
                                lon_min = lon - bufferlon,
                                lon_max = lon + bufferlon,
                                lat_min = lat - bufferlat,
                                lat_max = lat + bufferlat,
                                time_start = time - timedelta(hours=3),
                                time_end = time + duration)

o = OceanDrift()

o.add_reader(cmems)
o.seed_elements(lon=lon, lat=lat, number=5000, radius=1000, time=time)
o.run(duration=duration)
o.animation(fast=True, filename='cmems.gif')


#%%
# .. image:: /gallery/animations/cmems.gif
