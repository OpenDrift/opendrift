#!/usr/bin/env python
"""
Stokes drift vertical profiles
==============================
"""

import numpy as np
from opendrift.models.physics_methods import stokes_drift_profile_monochromatic, stokes_drift_profile_exponential, stokes_drift_profile_phillips, plot_stokes_profile

z = np.linspace(0, -10, 100)
su_surface = .3
sv_surface = .1
T = 8

#%%
# Plotting profiles of Stokes drift for significant wave heights of 2 meters
hs = 2
su_m, sv_m, ss_m = stokes_drift_profile_monochromatic(su_surface, sv_surface, hs, T, z)
su_e, sv_e, ss_e = stokes_drift_profile_exponential(su_surface, sv_surface, hs, T, z)
su_p, sv_p, ss_p = stokes_drift_profile_phillips(su_surface, sv_surface, hs, T, z)

profiles = [
    {'u': su_m, 'v': sv_m, 'z': z, 'kwargs': {'label': 'Stokes, Monochromatic'}},
    {'u': su_e, 'v': sv_e, 'z': z, 'kwargs': {'label': 'Stokes, Exponential'}},
    {'u': su_p, 'v': sv_p, 'z': z, 'kwargs': {'label': 'Stokes, Phillips'}}
    ]

plot_stokes_profile(profiles)
