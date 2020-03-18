#!/usr/bin/env python
"""
Plastic
==================================
"""

from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.plastdrift import PlastDrift

o = PlastDrift(loglevel=20)
o.list_configspec()  # to see available configuration options

# Arome atmospheric model
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
# Norkyst ocean model
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([reader_norkyst, reader_arome])
start_time = reader_arome.start_time
end_time = reader_arome.start_time + timedelta(hours=5)
end_time = reader_arome.end_time
time = [start_time, start_time + timedelta(hours=5)]

#%%
# Seeding some particles
lon = 4.6; lat = 60.0; # Outside Bergen
o.seed_elements(lon, lat, radius=50, number=3000, time=time)
o.run(end_time=end_time, time_step=1800, time_step_output=3600)

#%%
# Second run, without wind/Stokes drift
o2 = PlastDrift(loglevel=20)
o2.add_reader([reader_norkyst])
o2.seed_elements(lon, lat, radius=50, number=3000, time=time)
o2.run(end_time=end_time, time_step=1800, time_step_output=3600)

#%%
# Print and plot results
print(o)
o.animation(compare=o2, fast=True,
            legend=['Current + wind/Stokes drift', 'Current only'])
#o.animation(color='depth')
#o.plot_property('depth')

#%%
# .. image:: /gallery/animations/example_plast_0.gif

