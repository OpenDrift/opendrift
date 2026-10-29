#!/usr/bin/env python
"""
OceanDrift with NorKyst current and wind: THREDDS vs x3d
========================================================

Runs the same OceanDrift simulation (10 objects, starting now) with the NorKyst-800
surface current and wind from

1. the THREDDS OPeNDAP aggregate (standard), and
2. the x3d Zarr store on rgw.met.no, read through the x3d xarray backend
   (needs ``pip install x3d``),

and compares timings and trajectories.
"""

import time
from datetime import datetime, timedelta, timezone

import numpy as np
import xarray as xr
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

thredds_url = 'https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be'
x3d_url = ('x3d+https://rgw.met.no/5e4f32366fc74342866b4c67dbd4713d:'
           'norkyst800m-v3/norkyst800')

# The x3d store does not carry the projection variable of the THREDDS dataset
proj4 = ('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=70 +x_0=3369600 +y_0=1844800 '
         '+a=6378137 +b=6356752.3142 +units=m +no_defs')

start = datetime.now(timezone.utc).replace(tzinfo=None, minute=0, second=0, microsecond=0)

#%%
# Standard: NorKyst current and wind from THREDDS
t0 = time.perf_counter()
reader_thredds = reader_netCDF_CF_generic.Reader(thredds_url)
open_thredds = time.perf_counter() - t0

o_thredds = OceanDrift(loglevel=50)
o_thredds.add_reader(reader_thredds)
np.random.seed(0)  # identical initial positions in both runs
o_thredds.seed_elements(lon=4.5, lat=59.6, radius=100, number=10, time=start)

t0 = time.perf_counter()
ds_thredds = o_thredds.run(duration=timedelta(hours=24), time_step=1800,
                           time_step_output=3600)
run_thredds = time.perf_counter() - t0
print(f'THREDDS: open {open_thredds:.2f} s, run {run_thredds:.2f} s')

#%%
# Same data through the x3d xarray backend. The opened dataset is passed
# directly to the reader.
t0 = time.perf_counter()
reader_x3d = reader_netCDF_CF_generic.Reader(xr.open_dataset(x3d_url),
                                             proj4=proj4, name='norkyst_x3d')
open_x3d = time.perf_counter() - t0

o_x3d = OceanDrift(loglevel=50)
o_x3d.add_reader(reader_x3d)
np.random.seed(0)
o_x3d.seed_elements(lon=4.5, lat=59.6, radius=100, number=10, time=start)

t0 = time.perf_counter()
ds_x3d = o_x3d.run(duration=timedelta(hours=24), time_step=1800,
                   time_step_output=3600)
run_x3d = time.perf_counter() - t0
print(f'x3d:     open {open_x3d:.2f} s, run {run_x3d:.2f} s')

#%%
# Compare timing and trajectories. x3d stores the data as int16 (current 1 mm/s,
# wind 1 cm/s precision), so trajectories are close, but not identical.
print(f'Speedup (total): {(open_thredds + run_thredds) / (open_x3d + run_x3d):.1f}x')
print(f'Speedup (run):   {run_thredds / run_x3d:.1f}x')
dist = np.hypot(ds_thredds.lon.values - ds_x3d.lon.values,
                ds_thredds.lat.values - ds_x3d.lat.values) * 111e3
print(f'Max position difference: {np.nanmax(dist):.1f} m (approx.)')

#%%
# Plot both runs
o_thredds.plot(compare=o_x3d, legend=['THREDDS', 'x3d'], fast=True)
