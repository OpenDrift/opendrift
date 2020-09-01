#!/usr/bin/env python
"""
Ice berg (openberg) deterministic
==================================
"""

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openberg import OpenBerg

#%%
# Initialize model
steps = 60  # This is the number of forecast steps

o = OpenBerg()  # Basic drift model suitable for icebergs

#%%
# Preparing Readers
reader_current = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
   '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

reader_wind = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
   '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

o.add_reader([reader_current, reader_wind])

#%%
# Seeding elements
#
# Icebergs are moved with the ocean current as per Barker et al (2004),
# in addition to a fraction of the wind speed (wind_drift_factor).
# This factor depends on the properties of the elements.
# Default empirical values are:
# - Wind drift fraction: 0.018 (1.8 %) (Garret 1985)
# - Iceberg size: 	Keel dept = 60m
#					Waterline length = 90.5m


o.seed_elements(3.3, 60., radius=3000,
                time=reader_current.start_time,
                water_line_length=100,
                keel_depth=90, number=200)


print('Starting free run .../n')

print('Start time: ' + str(o.start_time))

#%%
# Running model

o.run(time_step=3600, steps=steps)

#%%
# Print and plot results
o.plot(fast=True)
o.animation(fast=True)


#%%
# .. image:: /gallery/animations/example_openberg_det_0.gif

print('############## Latitudes:', o.history['lat'])
print('############## Longitudes:',o.history['lon'])

