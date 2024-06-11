#!/usr/bin/env python
"""
Leeway
==================================
"""

from datetime import timedelta
import cmocean
import xarray as xr
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.leeway import Leeway

o = Leeway(loglevel=0)  # Set loglevel to 0 for debug information

# Atmospheric model for wind
#reader_arome = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Ocean model for current
#reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader(reader_norkyst)
o.add_reader(reader_arome)
o.set_config('environment:fallback:x_sea_water_velocity', 0)
o.set_config('environment:fallback:y_sea_water_velocity', 0)

#%%
# Seed leeway elements at defined position and time
object_type = 26  # 26 = Life-raft, no ballast
o.seed_elements(lon=4.5, lat=59.6, radius=100, number=1000,
                time=reader_arome.start_time, object_type=object_type)

#%%
# Running model
o.run(duration=timedelta(hours=48), time_step=900, time_step_output=3600)

#%%
# Print and plot results
print(o)

#%%
# Animation with current as background.
# Note that drift is also depending on wind, which is not shown.
o.animation(background=['x_sea_water_velocity', 'y_sea_water_velocity'],
             skip=5,  # show every 5th vector
             cmap=cmocean.cm.speed, vmin=0, vmax=.8, bgalpha=.7, land_color='#666666', fast=True)

#%%
# .. image:: /gallery/animations/example_leeway_0.gif

o.plot(fast=True)

#%%
# Plot density of stranded elements
d, dsub, dstr, lon, lat = o.get_density_array(pixelsize_m=3000)
strand_density = xr.DataArray(dstr[-1,:,:], coords={'lon_bin': lon[0:-1], 'lat_bin': lat[0:-1]})
o.plot(fast=True, background=strand_density.where(strand_density>0),
       vmin=0, vmax=20, cmap=cmocean.cm.dense, clabel='Density of stranded elements',
       show_elements=False, linewidth=0)
