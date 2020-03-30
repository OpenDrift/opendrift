#!/usr/bin/env python
"""
Diffusion
=============
"""

from datetime import datetime, timedelta

from opendrift.readers import reader_ECOM_building
from opendrift.models.oceandrift import OceanDrift


o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
reader_pcse = reader_ECOM_building.Reader('/home/arian/IC/arquivos_cdf/six_days.nc')

o.add_reader([reader_pcse])

		# Seeding some particles
lon = -44; lat = -24;

time = reader_pcse.start_time
o.seed_elements(lon, lat, radius=100, number=100, time=time)
#o.set_config('drift:current_uncertainty', 0)  # 0 is default
o.run(duration=timedelta(hours=100))

o.animation()
