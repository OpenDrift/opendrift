#!/usr/bin/env python
"""
Use a shapefile as landmask
===========================
"""

import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_shape
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

# This example works better using hourly input from Thredds than the daily data from test folder
reader_topaz = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/cmems/topaz6/dataset-topaz6-arc-15min-3km-be.ncml')

#%%
# Use shapes from Cartopy as landmask. These shapefiles are less acurate than those
# provided by the GSHHG dataset (available though the reader_global_landmask reader).
import cartopy.io.shapereader as shpreader
shpfilename = shpreader.natural_earth(resolution='110m',
                                    category='cultural',
                                    name='admin_0_countries')
reader_natural = reader_shape.Reader.from_shpfiles(shpfilename)

o.add_reader([reader_natural, reader_topaz])
o.set_config('general:use_auto_landmask', False)  # Disabling the automatic GSHHG landmask
o.set_config('general:coastline_action', 'stranding')

#%%
# Seeding elements on a grid
lons = np.linspace(12, 14.7, 30)
lats = np.linspace(67.5, 68.5, 30)
lons, lats = np.meshgrid(lons, lats)
lon = lons.ravel()
lat = lats.ravel()

time = reader_topaz.start_time
o.seed_elements(lon, lat, radius=0, number=30*30, time=time)

o.run(steps=48, time_step=3600)

#%%
# The custom landmask from shapefile, which is used for the simulation, is also shown in the map instead of the default GSHHG
o.plot()

o.animation()

#%%
# .. image:: /gallery/animations/example_shape_landmask_0.gif
