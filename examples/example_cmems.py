#!/usr/bin/env python

# This example runs a Leeway simulation, with current data downloaded
# from CMEMS, and wind data from an NCEP Thredds server.
# To run this example, you need a CMEMS account created at
# http://marine.copernicus.eu, and you need to install the motu_client.py
# available from https://github.com/clstoulouse/motu-client-python

import os
from datetime import datetime, timedelta
import numpy as np

from opendrift.models.leeway import Leeway
from opendrift.readers import reader_cmems
from opendrift.readers import reader_netCDF_CF_generic

try:
    f = open('cmems_account.txt')
except:
    raise ValueError('Please store your CMEMS username and password '
                     'in local file "cmems_account.txt", optionally '
                     'followed by a line with the path to motu_client.py')

username = f.readline().strip()
password = f.readline().strip()
try:
    motu_client = f.readline().strip()
except:
    motu_client = 'motu-client.py'
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
                                motu=motu_client,
                                lon_min = lon - bufferlon,
                                lon_max = lon + bufferlon,
                                lat_min = lat - bufferlat,
                                lat_max = lat + bufferlat,
                                time_start = time,
                                time_end = time + duration)

# Fetching wind data from NCEP
reader_ncep = reader_netCDF_CF_generic.Reader('http://oos.soest.hawaii.edu/thredds/dodsC/hioos/model/atm/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd')

o = Leeway()
o.add_reader([cmems, reader_ncep])
o.seed_elements(lon=lon, lat=lat, number=5000, time=time)
o.run(duration=duration, outfile='cmems.nc',
      time_step=600, time_step_output=3600)
o.animation()
