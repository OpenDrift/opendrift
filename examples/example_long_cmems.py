#!/usr/bin/env python
"""
Copernicus Marine Client (CMEMS)
================================

This example runs an OceanDrift simulation with current data from CMEMS.
To run this example, you need a CMEMS account created at https://marine.copernicus.eu
with username and password stored as environment variables
``COPERNICUSMARINE_SERVICE_USERNAME`` and ``COPERNICUSMARINE_SERVICE_PASSWORD``
or in a ``.netrc`` file with contents::

    machine copernicusmarine login <your username> password <your password>

This file must be stored in your home folder (and unreadable by others) or in the main OpenDrift folder

Alternatively, an Xarray dataset can be created explicitly with the copernicusmarine client, and provided to reader_netCDF_CF_generic:
https://opendrift.github.io/gallery/example_long_cmems_currents.html
"""

from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift

lon = 4.8; lat = 60  # Bergen, Norway
lon = -89; lat = 29.8  # New Orleans
lon = 107; lat = 10  # Ho Chi Minh
lon = 123; lat = -16.3  # Australia


o = OceanDrift()

o.add_readers_from_list(['cmems_mod_glo_phy_anfc_merged-uv_PT1H-i'])

o.seed_elements(lon=lon, lat=lat, number=5000, radius=1000, time=datetime.utcnow())
o.run(duration=timedelta(days=3))

o.animation(fast=True, clabel='Ocean current [m/s]',
            background=['x_sea_water_velocity', 'y_sea_water_velocity'])
