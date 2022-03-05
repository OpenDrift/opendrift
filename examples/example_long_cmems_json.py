#!/usr/bin/env python
"""
Initialising reader from JSON string
====================================
"""

from datetime import datetime, timedelta
from opendrift.models.leeway import Leeway

o = Leeway()

#%% 
# Adding a CMEMS reader as a JSON string, and NCEP winds from thredds URL.
# For CMEMS, username or password must be stored in a .netrc file under
# machine name "cmems", or provided explicitly instead of *null*
o.add_readers_from_list(['https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_anfc_merged-uv_PT1H-i',
    'https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd'])

o.seed_elements(time=datetime.utcnow(), lon=123, lat=-16.1,
                number=1000, radius=1000)

o.run(duration=timedelta(hours=72))

print(o)

o.animation(fast=False, skip=1, scale=10, background=
            ['x_sea_water_velocity', 'y_sea_water_velocity'])
