#!/usr/bin/env python
"""
Openberg - statistical mode
==============================
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_current_from_track
from opendrift.models.openberg import OpenBerg

#%%
# Create observation of iceberg track
obslon = [3.1, 3.123456]
obslat = [61.1, 61.132198]
obstime = [datetime(2015, 11, 16, 0), datetime(2015, 11, 16, 6)]

#%%
# Initialize model
steps = 60   # This is the number of forecast steps
o = OpenBerg(loglevel=30)  # Basic drift model suitable for icebergs

#%%
# Preparing Readers
reader_wind = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
   '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc',name='WIND')

reader_current = reader_current_from_track.Reader(
    obslon, obslat, obstime, wind_east=0, wind_north=0,
    windreader=reader_wind, wind_factor=0.018)

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
# 					NB! Iceberg size is irrelevant for current_reader with 1D z-profile

o.seed_elements(3.3, 61.3, radius=3000, number=500,
                time=reader_current.start_time)

#%%
# Run model
print('Starting free run .../n')

print('Start time: ' + str(o.start_time))

o.run(time_step=3600, steps=steps)

#%%
# Print and plot results
o.plot(fast=True)
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_openberg_stat_0.gif

