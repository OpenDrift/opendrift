#!/usr/bin/env python
"""
Ship drift
==================================
"""

from datetime import datetime
from opendrift.models.shipdrift import ShipDrift

o = ShipDrift(loglevel=20)

o.add_readers_from_list([
    'https://thredds.met.no/thredds/dodsC/cmems/topaz6/dataset-topaz6-arc-15min-3km-be.ncml',
    'https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc',
    'https://thredds.met.no/thredds/dodsC/cmems/mywavewam3km/dataset-wam-arctic-1hr3km-be.ncml'
    ])

#%%
# Seed ship elements at defined position and time
# Note: beam/length ratio is larger than allowed, but is then clipped internally
o.seed_elements(lon=5.0, lat=63.0, radius=1000, number=1000,
                time=datetime.utcnow(),
                length=80.0, beam=20.0, height=9.0, draft=4.0)

#%%
# Running model
o.run(steps=24, stop_on_error=True)

#%%
# Print and plot results
print(o)
o.plot(linecolor='orientation')
#o.animation(color='orientation', markersize=20, filename="orientation.gif", legend=['left','right'],colorbar=False,cmap='bwr')
o.animation(color='orientation', legend=['left','right'], markersize=20, colorbar=False, cmap='bwr')

#%%
# .. image:: /gallery/animations/example_shipdrift_0.gif
