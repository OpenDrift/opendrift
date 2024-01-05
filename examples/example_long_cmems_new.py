#!/usr/bin/env python
"""
Copernicus marine client
========================

This example runs an OceanDrift simulation
with current data from CMEMS
To run this example, you need a CMEMS account created at
https://marine.copernicus.eu
and the copernicus_marine_client installed from
https://pypi.org/project/copernicus-marine-client/
"""

from datetime import datetime, timedelta
import copernicus_marine_client as copernicusmarine
from opendrift.readers.reader_netCDF_CF_generic import Reader
from opendrift.models.oceandrift import OceanDrift

lon = 4.8; lat = 60  # Bergen, Norway
lon = -89; lat = 29.8  # New Orleans
lon = 107; lat = 10  # Ho Chi Minh
lon = 123; lat = -16.3  # Australia


o = OceanDrift()

#%%
# First get a Xarray dataset from copernicus_marine_client
ds = copernicusmarine.open_dataset(
    dataset_id='cmems_mod_glo_phy_anfc_merged-uv_PT1H-i',
    username='<your cmems username>', password='<your cmems password>')

#%%
# Then create an OpenDrift reader from this dataset
r = Reader(ds)
o.add_reader(r)

o.seed_elements(lon=lon, lat=lat, number=5000, radius=1000, time=datetime.utcnow())
o.run(duration=timedelta(days=3))

o.animation(fast=True, clabel='Ocean current [m/s]',
            background=['x_sea_water_velocity', 'y_sea_water_velocity'])

#%%
# .. image:: /gallery/animations/example_cmems_new_0.gif
