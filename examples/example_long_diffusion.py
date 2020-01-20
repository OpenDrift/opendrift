#!/usr/bin/env python
"""
Diffusion
=============
"""

from datetime import datetime, timedelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

lon = 4.5; lat = 60.0; # Outside Bergen

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Arome atmospheric model
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
# Norkyst ocean model
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
# Uncomment to use live data from thredds
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader([reader_norkyst, reader_arome])
time = reader_arome.start_time
o.seed_elements(lon, lat, radius=500, number=2000, time=time)
o.set_config('drift:current_uncertainty', 0)  # 0 is default
o.run(duration=timedelta(hours=24))

# Second run, identical, except for added diffusion

o2 = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
o2.add_reader([reader_norkyst, reader_arome])
o2.seed_elements(lon, lat, radius=500, number=2000, time=time)
o2.set_config('drift:current_uncertainty', .2) # Difference from first run
o2.run(duration=timedelta(hours=24))

# Comparing
o2.animation(filename='diffusion.gif', compare=o, legend=['0.2 m/s std for current components', 'No diffusion'], legend_loc='upper center')

#%%
# .. image:: /gallery/animations/diffusion.gif
