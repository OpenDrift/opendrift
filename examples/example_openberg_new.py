#!/usr/bin/env python
"""
Ice bergs (openberg)
====================
"""

from datetime import datetime, timedelta
from opendrift.models.openberg import OpenBerg

o = OpenBerg()

o.add_readers_from_list([
        'https://thredds.met.no/thredds/dodsC/cmems/topaz6/dataset-topaz6-arc-15min-3km-be.ncml',
        'https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd'])

o.seed_elements(lon=28, lat=77, time=datetime.now(), number=1000, radius=1000)

o.run(duration=timedelta(days=2))

o.plot(fast=True)
