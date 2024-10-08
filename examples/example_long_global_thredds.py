#!/usr/bin/env python
"""
Thredds (global)
==================================
"""

from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift

# Basic ocean drift module: current + 2% of wind
o = OceanDrift()

# Adding readers for global Thredds datasets:
# - Ocean forecast from global HYCOM model
# - Weather forecast from NOAA/NCEP
o.add_readers_from_list([
    'https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_uv3z/FMRC_ESPC-D-V02_uv3z_best.ncd',
    'https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd'])

o.seed_elements(lat=24, lon=-81, time=datetime.utcnow(),
                number=5000, radius=10000)

o.run(time_step=timedelta(minutes=30),
      duration=timedelta(days=5))
o.animation()
