#!/usr/bin/env python
"""
Diffusion
=============
"""

from datetime import datetime, timedelta

from opendrift.readers import reader_ECOM_building
from opendrift.models.oceandrift import OceanDrift


o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
reader_pcse = reader_ECOM_building.Reader(r'C:\Arian\arquivos_cdf\six_days.nc')
#reader_pcse= reader_netCDF_CF_generic.Reader(r'C:\Users\arian\opendrift\arquivos_netcdf\PCSE_300_u&v.nc')

o.add_reader([reader_pcse])

		# Seeding some particles
lon = -45.1; lat = -24;

time = reader_pcse.start_time
o.seed_elements(lon, lat, radius=500, number=2000, time=time)
o.set_config('drift:current_uncertainty', 0)  # 0 is default
o.run(duration=timedelta(hours=100))

#%%
# Second run, identical, except for added diffusion

o2 = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
o2.add_reader([reader_norkyst, reader_arome])
o2.seed_elements(lon, lat, radius=500, number=2000, time=time)
o2.set_config('drift:current_uncertainty', .2) # Difference from first run
o2.run(duration=timedelta(hours=100))

#%%
# Comparing
o2.animation(compare=o, legend=['0.2 m/s std for current components', 'No diffusion'], legend_loc='upper center')
