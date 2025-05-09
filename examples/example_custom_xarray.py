#!/usr/bin/env python
"""
Customising Xarray Dataset
==========================
"""

from datetime import datetime, timedelta
import xarray as xr
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Opening the currents netCDF file with xarray
ds = xr.open_dataset(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

#%%
# Creating and adding a landmask DataArray (variable) based on the u-current component
ds['landmask'] = ds.u.isel(time=0).isel(depth=0) * 0  # 0 (ocean) where current is finite
ds['landmask'] = ds.landmask.fillna(1)  # 1 (land) where current is NaN
ds['landmask'] = ds.landmask.assign_attrs(standard_name='land_binary_mask')  # Adding attribute standard_name so that this variable is recognised as landmask

#%%
# Creating an OpenDrift reader from this modified xarray dataset, and confirming that landmask is recognised
r = reader_netCDF_CF_generic.Reader(ds)
print(r)
o.add_reader(r)

o.set_config('general:use_auto_landmask', False)  # Disable the automatic GSHHG landmask, so that the custom landmask is used

#%%
# Seeding a particle and running simulation, and confirming that the custom landmask is used for stranding
o.seed_elements(lon=4.9, lat=60.0, time=r.start_time)
o.run(end_time=r.end_time)
o.plot(fast=True)
