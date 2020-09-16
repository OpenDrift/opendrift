#!/usr/bin/env python
"""
CMEMS
=============

This example runs an OceanDrift simulation
with current data downloaded from CMEMS
To run this example, you need a CMEMS account created at
https://marine.copernicus.eu, and you need to install the motuclient:
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
    print('Please store your CMEMS username and password '
          'in local file "cmems_account.txt"')
    quit()

cmems_user = f.readline().strip()
cmems_password = f.readline().strip()
f.close()

#%%
# Seed information
lon = 4.8; lat = 60  # Bergen, Norway
lon = -89; lat = 29.8  # New Orleans
lon = 107; lat = 10  # Ho Chi Minh
lon = 123; lat = -16.3  # Australia


time = datetime.utcnow()
duration = timedelta(days=3)
bufferlat = duration.total_seconds()/111000
bufferlon = bufferlat*np.cos(lat*np.pi/180)

#%% Fetching current data from CMEMS
# By intialisation, only an XML-file with the contents is downloaded.
# The netCDF file with data is downloaded later, after the simulation is started,
# when necessary coverage (time/space) is known ( <reader>.prepare() ) 
cmems = reader_cmems.Reader(
    dataset='global-analysis-forecast-phy-001-024-hourly-merged-uv',
    variable_mapping={  # Overriding the mapping in reader_cmems.py
        'utotal': 'x_sea_water_velocity',
        'vtotal': 'y_sea_water_velocity',
        'utide': 'sea_ice_x_velocity',   # Fake mapping, as standard_name
        'vtide': 'sea_ice_y_velocity'},  # must be unique for a reader
    cmems_user=cmems_user, cmems_password=cmems_password)

o = OceanDrift()

o.add_reader(cmems)
o.seed_elements(lon=lon, lat=lat, number=5000, radius=1000, time=time)
o.run(duration=duration)

# Although total current (SMOC) has been used as forcing, we plot the tidal current as background field, disguised as ice velocity
o.animation(fast=True, clabel='Tidal current [m/s]', skip=1, scale=20, 
            background=[ 'sea_ice_x_velocity', 'sea_ice_y_velocity'])

#%%
# .. image:: /gallery/animations/example_cmems_0.gif
