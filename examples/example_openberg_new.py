#!/usr/bin/env python
"""
Icebergs (openberg)
====================
"""

from opendrift.models.openberg import OpenBerg
from datetime import datetime,timedelta
 
o = OpenBerg(add_stokes_drift = True,
             wave_rad = True,
             grounding = False,
             vertical_profile = False,
             melting = False,
             choose_melting  = {"wave": True, "lateral": True, "basal": True})

o.add_readers_from_list([
        'https://thredds.met.no/thredds/dodsC/cmems/topaz6/dataset-topaz6-arc-15min-3km-be.ncml',
        'https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd'])


o.seed_elements(lon= -57,lat= 69,time=datetime.now(), number=1000, radius=1000)
o.run(duration=timedelta(days=2)) 
o.plot(fast=True)
