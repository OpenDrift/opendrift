"""
Combining readers using operators
==================================

It is possible to combine readers using operators, to create e.g. a mean reader from different sources, or adding a constant force term to another sources.

.. seealso::

    :py:mod:`opendrift.readers.operators`.

"""

from datetime import timedelta
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
# First a regular simulation, without the added wind:
lw.add_reader([reader_arome, reader_norkyst])
object_type = 26  # 26 = Life-raft, no ballast
lw.seed_elements(lon=4.5, lat=59.6, radius=100, number=1000,
                 time=reader_arome.start_time, object_type=object_type)
lw.run(duration=timedelta(hours=48), time_step=900, time_step_output=3600)
lw.plot(fast=True)

#%%
# Then a simulation where we add a constant x_wind component to cause stranding.
reader_onshore = reader_constant.Reader({ 'x_wind' : 10., 'y_wind': 0 })

#%%
# If we just add all readers OpenDrift will read the wind from the first available reader, so we combine the wind reader and add them in order so
# that the combined reader is used first, then the actual sources are used next.
r0 = reader_onshore + reader_arome
print(r0)

lw2 = Leeway(loglevel=20)  # Set loglevel to 0 for debug information
lw2.add_reader([r0, reader_arome, reader_norkyst])
#%%
# Seed leeway elements at defined position and time
object_type = 26  # 26 = Life-raft, no ballast
lw2.seed_elements(lon=4.5, lat=59.6, radius=100, number=1000,
                 time=reader_arome.start_time, object_type=object_type)

#%%
# Running model
lw2.run(duration=timedelta(hours=48), time_step=900, time_step_output=3600)

#%%
# Print and plot results
print(lw2)
lw2.plot(fast=True)

#%%
# The first simulation compared to the second with extra wind:
lw.plot(fast=True, compare=lw2)

