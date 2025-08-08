#!/usr/bin/env python
"""
Caspian / lake
===============
"""

from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_shape
import cartopy

#%%
# The default GSHHG landmask reader classifies lakes
# (e.g. Caspian Sea) as land and not water.
# For simulations in lakes, we can use a custom reader for lakes only.
level = 'h'  # using high resolution coastline
reader_lakes = reader_shape.Reader.from_shpfiles(
    f'{cartopy.config['data_dir']}/shapefiles/gshhs/{level}/GSHHS_{level}_L2.shp',
    invert=True)  # Inverting since inside of polygons is water (lake) and not land.

#%%
# Disabling the default landmask reader, and using the above reader instead
o = OceanDrift(loglevel=20)
o.add_reader(reader_lakes)
o.set_config('general:use_auto_landmask', False)  # To use custom landmask instead
o.set_config('environment:constant:x_sea_water_velocity', 1)
o.set_config('drift:horizontal_diffusivity', 100)
o.seed_elements(lon=48.819, lat=44.959, radius=5000, number=100, time=datetime.now())

o.run(steps=10)
o.plot(fast=False, buffer=2)
