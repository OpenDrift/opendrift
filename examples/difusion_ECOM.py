#!/usr/bin/env python
"""
Diffusion_ECOM_EC1 (Verão 2014)
=============
"""

from datetime import datetime, timedelta

from opendrift.readers import reader_ECOM_building
from opendrift.models.oceandrift import OceanDrift


o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
reader_pcse = reader_ECOM_building.Reader('/home/arian/IC/arquivos_cdf/six_days.nc')

o.add_reader([reader_pcse])

		# Seeding some particles
lon = -45.6; lat = -23.9;

time = reader_pcse.start_time
o.seed_elements(lon, lat, radius=100, number=200,time=time)
#o.set_config('drift:current_uncertainty', 0)  # 0 is default
#o.run(duration=timedelta(hours=120))
o.run(steps=66*2, time_step=1800)
#o.animation()

o2 = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
o2.add_reader([reader_pcse])
o2.seed_elements(lon, lat, radius=100, number=200, time=time)
o2.set_config('drift:current_uncertainty', .09) # Difference from first run
#o2.run(duration=timedelta(hours=120))
o2.run(steps=66*2, time_step=1800)
#%%
# Comparing
o2.plot(fast=True)
o2.animation(compare=o, legend=['0.2 m/s std for current components', 'No diffusion'], legend_loc='upper center', filename ='Diffusion_ECOM_2014.gif' )
