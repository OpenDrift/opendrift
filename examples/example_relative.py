#!/usr/bin/env python
"""
Relative and absolute wind
==================================
"""

from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader([reader_norkyst, reader_arome])

#%%
# Seeding some particles
lon = 4.2; lat = 60.0; # Outside Bergen
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]

#%%
# Using windspeed relative to moving ocean (current)
o.set_config('drift:relative_wind',  False)
o.set_config('drift:vertical_mixing', False)

#%%
# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=50, number=5000, time=time)

o.run(steps=48*2, time_step=1800, time_step_output=3600*2)

#%%
# Second run, for comparison
o2 = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o2.add_reader([reader_norkyst, reader_arome])
o2.set_config('drift:relative_wind',  True)
o2.set_config('drift:vertical_mixing', False)
o2.seed_elements(lon, lat, radius=50, number=5000, time=time)
o2.run(steps=48*2, time_step=1800, time_step_output=3600*2)


#%%
# Animate and compare the two runs
o.animation(compare=o2, legend=['Absolute wind', 'Relative wind'])

#%%
# .. image:: /gallery/animations/example_relative_0.gif

o.plot(compare=o2, legend=['Absolute wind', 'Relative wind'], fast=True)
