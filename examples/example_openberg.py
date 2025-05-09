#!/usr/bin/env python
"""
Icebergs (openberg)
====================
"""

from opendrift.models.openberg import OpenBerg
from datetime import datetime,timedelta
 
o = OpenBerg()

# The user can overwrite the default setup using set_config method
o.set_config('drift:vertical_profile', False) # use surface currents for this test

o.add_readers_from_list([
        'https://thredds.met.no/thredds/dodsC/cmems/topaz6/dataset-topaz6-arc-15min-3km-be.ncml',
        'https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd'])


o.seed_elements(lon=-56, lat=72, time=datetime.now(),
                number=100, radius=500,
                sail=10,draft=50,length=90,width=40)

o.run(duration=timedelta(days=3))
o.plot(fast=True)
o.plot_property('draft')