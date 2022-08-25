"""
Combining readers using operators
==================================

It is possible to combine readers using operators, to create e.g. a mean reader from different sources, or adding a constant force term to another sources.
"""

from datetime import timedelta
import cmocean
import xarray as xr
from opendrift.readers import reader_netCDF_CF_generic, reader_constant
from opendrift.models.leeway import Leeway

lw = Leeway(loglevel=20)  # Set loglevel to 0 for debug information

# Atmospheric model for wind
#reader_arome = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')
reader_arome = reader_netCDF_CF_generic.Reader(lw.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Ocean model for current
#reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')
reader_norkyst = reader_netCDF_CF_generic.Reader(lw.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

#%%
# We add a constant x_wind component to cause stranding.
reader_onshore = reader_constant.Reader({ 'x_wind' : 10., 'y_wind': 0 })

#%%
# If we just add all readers OpenDrift will read the wind from the first available reader, so we combine the wind reader and add them in order so
# that the combined reader is used first, then the actual sources are used next.
r0 = reader_onshore + reader_arome
print(r0)

lw.add_reader([r0, reader_arome, reader_norkyst])
#%%
# Seed leeway elements at defined position and time
object_type = 26  # 26 = Life-raft, no ballast
lw.seed_elements(lon=4.5, lat=59.6, radius=100, number=1000,
                 time=reader_arome.start_time, object_type=object_type)

#%%
# Running model
lw.run(duration=timedelta(hours=48), time_step=900, time_step_output=3600)

#%%
# Print and plot results
print(lw)

#%%
# Animation with current as background.
# Note that drift is also depending on wind, which is not shown.
# lw.animation(background=['x_sea_water_velocity', 'y_sea_water_velocity'],
#              skip=5,  # show every 5th vector
#              cmap=cmocean.cm.speed, vmin=0, vmax=.8, bgalpha=.7, land_color='#666666', fast=True)

#%%
# .. image:: /gallery/animations/example_leeway_0.gif

lw.plot(fast=True)

