#!/usr/bin/env python
"""
Reader boundary
==================================

Seeding elements around the border of a ocean model domain (NorKyst800)
to demonstrate autmatic transition back and forth with another model
covering a larger domain (Nordic).
"""

import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

# Nordic4
reader_nordic4 = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/nordic4km/zdepths1h/aggregate_be')

o.add_reader([reader_norkyst, reader_nordic4])
o.fallback_values['land_binary_mask'] = 0

#%%
# Seeding some particles
lons = np.linspace(10.2, 12.2, 50)
lats = np.linspace(69.8, 70.8, 50)
lons, lats = np.meshgrid(lons, lats)

#%%
# Seed oil elements at defined position and time
o.seed_elements(lons, lats, radius=0, number=2500,
                time=reader_norkyst.start_time)

#%%
# Running model (until end of driver data)
o.run(steps=16*4, time_step=900, time_step_output=1800)

#%%
# Print and plot results
print(o)
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_reader_boundary_0.gif

o.plot(fast=True)
