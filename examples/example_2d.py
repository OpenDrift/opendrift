#!/usr/bin/env python
"""
2D simulation profile
=====================
"""

from datetime import datetime, timedelta
import numpy as np
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Disable any 3D motion
o.disable_vertical_motion()

#%%
# Define some constant current, wind and Stokes drift
o.set_config('environment:fallback:x_wind', 7)
o.set_config('environment:fallback:x_sea_water_velocity', .1)
o.set_config('environment:fallback:sea_surface_wave_stokes_drift_x_velocity', .2)
o.set_config('environment:fallback:sea_surface_wave_significant_height', 2)
o.set_config('environment:fallback:sea_surface_wave_period_at_variance_spectral_density_maximum', 8)

#%%
# Seed elements between surface and 5m depth
time = datetime.utcnow()
z = -np.linspace(0, 5, 50)
o.seed_elements(lon=4.5, lat=60.0, z=z, radius=0, number=len(z), time=time)

#%%
# Running model for 6 hours
o.run(duration=timedelta(hours=6), time_step=600)

#%%
# To check that z is really kept constant for all particles
o.plot_property('z')

#%%
# Vertical profile of advection. Note the decaying importance of Stokes drift, and the additional windage of the element at surface
o.animation_profile()

#%%
# .. image:: /gallery/animations/example_2d_0.gif
