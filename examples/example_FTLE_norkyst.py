#!/usr/bin/env python
"""
LCS Norkyst
==================================
"""

from datetime import datetime, timedelta
import numpy as np
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

# Norkyst ocean model
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

#o.add_reader([reader_norkyst, reader_arome])
o.add_reader([reader_norkyst])


#%%
# Calculating attracting/backwards FTLE at 20 hours
lcs = o.calculate_ftle(
    time=reader_norkyst.start_time + timedelta(hours=24),
    time_step=timedelta(minutes=15),
    duration=timedelta(hours=3), delta=800,
    RLCS=False)

#%%
# Simulation from beginning and up to 30 hours (time of LCS)
o.reset()
lons = np.linspace(3.2, 5.0, 100)
lats = np.linspace(59.8, 61, 100)
lons, lats = np.meshgrid(lons, lats)
lons = lons.ravel()
lats = lats.ravel()
o.seed_elements(lons, lats, radius=0, number=10000,
                time=reader_norkyst.start_time)

o.run(end_time=reader_norkyst.start_time+timedelta(hours=24),
      time_step=timedelta(minutes=30))

o.plot(lcs=lcs, vmin=1e-7, vmax=1e-4, cmap='Reds', markersize=1, colorbar=True, show_particles=True, show_initial=False, linewidth=0)
