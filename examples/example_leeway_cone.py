#!/usr/bin/env python

from datetime import timedelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.leeway import Leeway

lw = Leeway(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(lw.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(lw.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

lw.add_reader([reader_norkyst, reader_arome])
lw.fallback_values['x_sea_water_velocity'] = 0
lw.fallback_values['y_sea_water_velocity'] = 0
lw.fallback_values['x_wind'] = 0
lw.fallback_values['y_wind'] = 0

# Intermediate map resolution is sufficient for large scale
lw.set_config('general:basemap_resolution', 'i')

# Seed elements along cone, e.g. ship track with
# increasing uncertainty in position
lon = [3.6, 5.1]; lat = [61., 59.6];
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]

objType = 26  # 26 = Life-raft, no ballast
lw.seed_elements(lon, lat, radius=[1000, 10000], number=5000,
                 time=time, objectType=objType)

# Running model
lw.run(steps=66*4, time_step=900)
print(lw)

# Print and plot results
print(lw)
lw.plot()
lw.animation()
