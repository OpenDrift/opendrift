#!/usr/bin/env python
"""
Trajan demo
============

From OpenDrift 2.0, analysis and plotting of results from OpenDrift simulations
will be handled by a new, standalone package: Trajan
https://github.com/OpenDrift/trajan

This example creates a test dataset, and demonstrates its anlysis using Trajan
"""

import os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import xarray as xr
import trajan as ta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

#%%
# Create test dataset with OpenDrift

o = OpenOil(loglevel=20)

# Add forcing
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
o.add_reader([reader_norkyst, reader_arome])

# Seeding some particles
o.seed_elements(lon=4.4, lat=60.1, radius=1000, number=1000,
                time=reader_arome.start_time)

# Running model
o.run(end_time=reader_norkyst.end_time, outfile='openoil.nc')

#%%
# Demonstrating analysis and visualisation of the output dataset, independently of OpenDrift code

if not os.path.exists('openoil.nc'):
    raise ValueError('Please run create_test_dataset.py first')

#%%
# Importing a trajectory dataset from a simulation with OpenDrift
# decode_coords is needed so that lon and lat are not interpreted as coordinate variables
d = xr.open_dataset('openoil.nc', decode_coords=False)
# Requirement that status>=0 is needed since non-valid points are not masked in OpenDrift output
d = d.where(d.status>=0)  # only active particles


#%%
# Displaying a basic plot of trajectories
d.traj.plot(land='mask')
ax = plt.gca()
ax.set_title('Adding custom title')
plt.show()

#%%
# Demonstrating how the Xarray Dataset can be modified, allowing for
# more flexibility than can be provided through the plotting method of OpenDrift

#%%
# Extracting only the first 100 elements, and every 4th output time steps:
dsub = d.isel(trajectory=range(0, 100), time=range(0, len(d.time), 4))
dsub.traj.plot(land='h')
plt.show()

#%%
# With several plots on the same figure
d.traj.plot(color='red', alpha=0.01, land='mask')  # Plotting individual trajectories in red
dmean = d.mean('trajectory', skipna=True) # Overlaying a "mean" trajectory in black
dmean.traj.plot(color='k', linewidth=5)
# Showing the a sub-period of the mean trajectory in yellow
dmean.sel(time=slice('2015-11-17', '2015-11-17 12')).traj.plot(color='yellow', linewidth=5)
plt.tight_layout()
plt.show()

#%%
# Cleaning up
os.remove('openoil.nc')
