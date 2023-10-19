#!/usr/bin/env python
"""
No stranding
=============
"""

import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o.set_config('drift:max_speed', 3)

# This example works better using hourly input from Thredds than the daily data from test folder
reader_current = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/cmems/topaz6/dataset-topaz6-arc-15min-3km-be.ncml')
#reader_current = reader_ROMS_native.Reader(o.test_data_folder() +
#    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

o.add_reader([reader_current])
o.set_config('general:coastline_action', 'previous')

#%%
# Seeding elements on a grid
lons = np.linspace(12, 14.5, 30)
lats = np.linspace(67.5, 68.5, 30)
lons, lats = np.meshgrid(lons, lats)
lon = lons.ravel()
lat = lats.ravel()

time = reader_current.start_time
o.seed_elements(lon, lat, radius=0, number=30*30, time=time)

o.run(steps=24*2, time_step=3600)

#%%
# Print and plot results
print(o)
o.plot()
#o.plot(background=['x_sea_water_velocity', 'y_sea_water_velocity'])
o.animation()

#%%
# .. image:: /gallery/animations/example_coastline_0.gif
