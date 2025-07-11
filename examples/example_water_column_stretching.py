#!/usr/bin/env python
"""
Water column stretching
========================
"""

from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from opendrift.readers import reader_oscillating
from opendrift.models.oceandrift import OceanDrift


#%%
# In OpenDrift, the vertical position of elements ("z") is defined as relative to actual surface,
# and not an absolute reference level (e.g. mean sea surface height).
# Thus if sea surface elevation changes with time (e.g. tides),
# we need to add a "correction / perturbation" to z, otherwise elements at/near seafloor
# will be lifted if surface elevation increases and z (relative to surface) remains unchanged.
# This correction is presently only implemented for OceanDrift, and must be switched on
# with config setting "drift:water:column_stretching"

# To illustrate, we add a reader with oscillating sea surface elevation (tidal)
# with amplitude of 1m and peroid of 6 hours
time = datetime.now()
reader_tidal = reader_oscillating.Reader('sea_surface_height', amplitude=-1,
                                         period=timedelta(hours=6), zero_time=time)

#%%
# First an illustration withouth this correction.
o = OceanDrift(loglevel=20)
o.add_reader(reader_tidal)
o.set_config('drift:water_column_stretching', False)
o.set_config('environment:constant:sea_floor_depth_below_sea_level', 10)
z = np.arange(0, -11, -1)  # Seeding one particle every meter from surface to 10m depth
o.seed_elements(lon=0, lat=0, time=time, z=z, number=11)
o.run(duration=timedelta(hours=24), time_step=1800)
o.result.z.plot.line(x='time', add_legend=False)
plt.show()

#%%
# We see that the particles remain at their initial depths (since we have no vertical advection or mixing),
# except for the element starting at seafloor, which is lifted up when sea level rises,
# since the config setting `drift:seafloor_action` is `lift_to_seafloor` by default.
# This lifting is in this case unphysical.

#%%
# We then make a simulation wih correction for the stretching/contraction of the water column.
o = OceanDrift(loglevel=20)
o.add_reader(reader_tidal)
o.set_config('drift:water_column_stretching', True)
o.set_config('environment:constant:sea_floor_depth_below_sea_level', 10)
o.seed_elements(lon=0, lat=0, time=time, z=z, number=11)
o.run(duration=timedelta(hours=24), time_step=1800)
o.result.z.plot.line(x='time', add_legend=False)
plt.show()

#%%
# Here we see that element depth (z, relative to surface) is changed so that
# elements at surface and seafloor remain at resp surface (z=0) and
# seafloor (z = sea_floor_depth + sea_surface_elevation)
