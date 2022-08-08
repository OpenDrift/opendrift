#!/usr/bin/env python
"""
CMEMS
=============

This example runs an OceanDrift simulation
with current data from CMEMS
To run this example, you need a CMEMS account created at
https://marine.copernicus.eu
with username and password stored in a ``.netrc`` file with contents::

    machine nrt.cmems-du.eu login <your username> password <your password>

This file must be stored in your home folder and made unreadable by others with ``chmod 700 .netrc``
Additionally, a file ``.dodsrc`` must also be stored in your home folder, containing the following line::

    HTTP.NETRC=<path to your .netrc file>

"""

from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift

lon = 4.8; lat = 60  # Bergen, Norway
lon = -89; lat = 29.8  # New Orleans
lon = 107; lat = 10  # Ho Chi Minh
lon = 123; lat = -16.3  # Australia


o = OceanDrift()

o.add_readers_from_list(['https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_anfc_merged-uv_PT1H-i'])

o.seed_elements(lon=lon, lat=lat, number=5000, radius=1000, time=datetime.utcnow())
o.run(duration=timedelta(days=3))

o.animation(fast=True, clabel='Ocean current [m/s]',
            background=['x_sea_water_velocity', 'y_sea_water_velocity'])

#%%
# .. image:: /gallery/animations/example_cmems_0.gif
