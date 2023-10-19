#!/usr/bin/env python
"""
Horizontal diffusion
=====================
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

lon = 4.5; lat = 60.0; # Outside Bergen

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Adding readers

# Arome atmospheric model
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
# Norkyst ocean model
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
# Uncomment to use live data from thredds
#reader_arome = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader([reader_norkyst, reader_arome])

#%%
# First run, with no horizontal diffusion
o.set_config('drift:current_uncertainty', 0)
o.set_config('drift:wind_uncertainty', 0)
time = reader_arome.start_time
o.seed_elements(lon, lat, radius=500, number=2000, time=time)
o.run(duration=timedelta(hours=24))

#%%
# Second run, identical, except for added diffusion
o2 = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o2.add_reader([reader_norkyst, reader_arome])
#o2.set_config('drift:current_uncertainty', .2) # Difference from first run
#o2.set_config('drift:wind_uncertainty', 1)     # Difference from first run
o2.set_config('drift:horizontal_diffusivity', 10)     # Difference from first run
o2.seed_elements(lon, lat, radius=500, number=2000, time=time)
o2.run(duration=timedelta(hours=24))

#%%
# Third run, identical, except for diffusion and shorter timestep
o3 = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o3.add_reader([reader_norkyst, reader_arome])
#o3.set_config('drift:current_uncertainty', .2) # Difference from first run
#o3.set_config('drift:wind_uncertainty', 1)     # Difference from first run
o3.set_config('drift:horizontal_diffusivity', 10)     # Difference from first run
o3.seed_elements(lon, lat, radius=500, number=2000, time=time)
o3.run(duration=timedelta(hours=24), time_step=300, time_step_output=3600)

#%%
# Comparing
o2.animation(compare=[o3, o], legend=['Diffusion, timstep 3600s', 'Diffusion, timestep 300s', 'No diffusion'],
             legend_loc='upper center', fast=True)

#%%
# .. image:: /gallery/animations/example_horizontal_diffusion_0.gif
