#!/usr/bin/env python
"""
Use a shapefile as landmask
===========================
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_shape
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o.max_speed = 3

# This example works better using hourly input from Thredds than the daily data from test folder
reader_nordic = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/nordic4km/zdepths1h/aggregate_be')

#%%
# Use shapes from Cartopy for tests. These shapefiles are less acurate than those
# provided by the GSHHS dataset (available though the reader_global_landmask reader).
import cartopy.io.shapereader as shpreader
shpfilename = shpreader.natural_earth(resolution='110m',
                                    category='cultural',
                                    name='admin_0_countries')
reader_natural = reader_shape.Reader.from_shpfiles(shpfilename)
#reader_nordic = reader_ROMS_native.Reader(o.test_data_folder() +
#    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

o.add_reader([reader_natural, reader_nordic])
o.set_config('general:use_auto_landmask', False)
o.set_config('general:coastline_action', 'stranding')

#%%
# Seeding elements on a grid
lons = np.linspace(12, 14.5, 30)
lats = np.linspace(67.5, 68.5, 30)
lons, lats = np.meshgrid(lons, lats)
lon = lons.ravel()
lat = lats.ravel()

time = reader_nordic.start_time
o.seed_elements(lon, lat, radius=0, number=30*30, time=time)

o.run(steps=48*2, time_step=3600)

#%%
# Print and plot results
print(o)
ax, _ = o.plot(hide_landmask = True, show = False)

#%%
# Show shapes
ax.add_geometries(reader_natural.polys, ccrs.PlateCarree(), facecolor=cfeature.COLORS['land'], edgecolor='black')
plt.show()

o.animation()
#%%
# .. image:: /gallery/animations/example_shape_landmask_0.gif
