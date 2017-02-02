#!/usr/bin/env python

import os
import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native

from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
o.max_speed = 3

# This example works better using hourly input from Thredds than the daily data from test folder
#reader_nordic = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/fou-hi/nordic4km-1h/Nordic-4km_SURF_1h_avg_00.nc')
reader_nordic = reader_ROMS_native.Reader(o.test_data_folder() +
    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

# We do not want landmask from ocean model
reader_nordic.variables = [v for v in reader_nordic.variables
                           if v != 'land_binary_mask']

# linearND interpolates data towards coast; uncomment for nearest interpolation
reader_nordic.interpolation='linearND'

print reader_nordic

o.add_reader([reader_nordic])
o.set_config('general:basemap_resolution', 'h')
o.set_config('general:coastline_action', 'previous')

# Seeding elements on a grid
lons = np.linspace(12, 14.5, 30)
lats = np.linspace(67.5, 68.5, 30)
lons, lats = np.meshgrid(lons, lats)
lon = lons.ravel()
lat = lats.ravel()

time = reader_nordic.start_time
o.seed_elements(lon, lat, radius=0, number=30*30, time=time)

o.run(steps=24*2, time_step=3600)

# Print and plot results
print o
o.plot()
#o.plot(background=['x_sea_water_velocity', 'y_sea_water_velocity'])
o.animation()
