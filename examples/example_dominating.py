#!/usr/bin/env python
"""
Show dominating source
======================
"""

import os
from datetime import datetime
import xarray as xr
import opendrift
from opendrift.models.oceandrift import OceanDrift

of = 'test.nc'

#%%
# Seed elements at 5 different locations/longitudes
lons = [4, 4.2, 4.3, 4.32, 4.6]
t = datetime.now()

o = OceanDrift(loglevel=20)
o.set_config('environment:constant:y_sea_water_velocity', .1)

for i, lon in enumerate(lons):
    o.seed_elements(lon=lon, lat=60, radius=3000, number=2000, time=t, origin_marker_name='Lon %f' % lon)
o.run(steps=15, outfile=of)

#%%
# Calculate spatial density of elements at 1500m grid spacing
d = o.get_histogram(pixelsize_m=1500, weights=None)
dom = d.argmax(dim='origin_marker', skipna=True)
dom = dom.where(d.sum(dim='origin_marker')>0)
dom.name = 'Dominating source'

#%%
# Show which of the 5 sources are dominating within each grid cell
o.animation(background=dom, show_elements=False, bgalpha=1,
             legend=list(o.origin_marker.values()), colorbar=False, vmin=0, vmax=4)


#%%
# .. image:: /gallery/animations/example_dominating_0.gif

os.remove(of)
