#!/usr/bin/env python
"""
Compare
=============

Comparing two oil drift simulation runs, with and without wind
"""

from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

# Arome atmospheric model
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
# Norkyst ocean model
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([reader_norkyst, reader_arome])

# Seeding some particles
lon = 4.5; lat = 60.0; # Outside Bergen
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]
o.seed_elements(lon, lat, radius=50, number=5000, time=time,
                oil_type='GENERIC HEAVY CRUDE',
                wind_drift_factor=0.03) # 3% wind drift

# Running model
o.run(steps=66, time_step=1800, time_step_output=3600)

# Second run, for comparison
o2 = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information
o2.add_reader([reader_norkyst, reader_arome])
o2.seed_elements(lon, lat, radius=50, number=5000, time=time,
                 oil_type='GENERIC HEAVY CRUDE',
                 wind_drift_factor=0.0) # No wind drift
o2.run(steps=66, time_step=1800, time_step_output=3600)

#%%
# Animate and compare the two runs.
# We see that there is much more stranding of oil when wind is considered.
o.animation(fast=True, compare=o2,
            legend=['Current + 3 % wind drift', 'Current only'])

#%%
# .. image:: /gallery/animations/example_compare_0.gif

o.plot(fast=True, compare=o2,
       legend=['Current + 3 % wind drift', 'Current only'])
