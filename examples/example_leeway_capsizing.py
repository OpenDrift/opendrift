#!/usr/bin/env python
"""
Leeway capsizing
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

#%%
# Activating capsizing for high winds, with probability per hour given by
# p(windspeed) = 0.5 + 0.5*tanh((windspeed-wind_threshold)/sigma)
o.set_config('processes:capsizing', True)
o.set_config('capsizing:wind_threshold', 30)
o.set_config('capsizing:wind_threshold_sigma', 5)
o.set_config('capsizing:leeway_fraction', 0.4)  # Reducing leeway coefficients to 40% of original after capsize

o.plot_capsize_probability()

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
# Animation illustrating effect of capsizing
from matplotlib.colors import ListedColormap
o.animation(color='capsized', cmap=ListedColormap(['black','red']), fast=True)

#%%
# .. image:: /gallery/animations/example_leeway_capsizing_0.gif

#%%
# Reverse run, also with probability of (un)-capsizing
o = Leeway()
o.add_reader(reader_norkyst)
o.add_reader(reader_arome)
o.set_config('processes:capsizing', True)
o.seed_elements(lon=4.4, lat=61.0, radius=100, number=1000,
                capsized=1,  # now we seed all objects as already capsized
                time=reader_arome.end_time, object_type=object_type)
o.run(time_step=-900, duration=timedelta(hours=48), time_step_output=timedelta(hours=1))
o.animation(color='capsized', cmap=ListedColormap(['black','red']), fast=True)

#%%
# .. image:: /gallery/animations/example_leeway_capsizing_1.gif
