#!/usr/bin/env python
"""
Stokes drift vertical profiles
==============================
"""

import numpy as np
from opendrift.models.physics_methods import stokes_drift_profile_breivik, plot_stokes_profile

z = np.linspace(0, -10, 100)
su_surface = .3
sv_surface = .1
T = 8

#%%
# Plotting profiles of Stokes drift for significant wave heights of 1 and 2 meters

hs = 1
su1, sv1, ss1 = stokes_drift_profile_breivik(su_surface, sv_surface, hs, T, z)
hs = 2
su2, sv2, ss2 = stokes_drift_profile_breivik(su_surface, sv_surface, hs, T, z)

profiles = [
    {'u': su1, 'v': sv1, 'z': z, 'kwargs': {'label': 'Stokes, Hs = 1m'}},
    {'u': su2, 'v': sv2, 'z': z, 'kwargs': {'label': 'Stokes, Hs = 2m'}}
    ]

plot_stokes_profile(profiles)
