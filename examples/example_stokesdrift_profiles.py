#!/usr/bin/env python
"""
Stokes drift vertical profiles
==============================
"""

import numpy as np
from opendrift.models.physics_methods import stokes_drift_profile_monochromatic, stokes_drift_profile_exponential, \
                                             stokes_drift_profile_phillips, stokes_drift_profile_windsea_swell, plot_stokes_profile

z = np.linspace(0, -10, 1000)

#%%
# Surface stokes drift and wave conditions
su_surface = .3
sv_surface = .1
Tp = 8
hs = 2
swell_dir_to = 0  # northwards
swell_mean_period = 10  # seconds
swell_height = 1.5
wind_sea_height = np.sqrt(hs**2 - swell_height**2)
wind_sea_mean_dir_to = np.degrees(np.arctan2(su_surface, sv_surface)) + 30
wind_sea_mean_period = 2  # seconds

#%%
# Plotting profiles of Stokes drift for significant wave heights of 2 meters
su_m, sv_m, ss_m = stokes_drift_profile_monochromatic(su_surface, sv_surface, hs, Tp, z)
su_e, sv_e, ss_e = stokes_drift_profile_exponential(su_surface, sv_surface, hs, Tp, z)
su_p, sv_p, ss_p = stokes_drift_profile_phillips(su_surface, sv_surface, hs, Tp, z)
su_c, sv_c, ss_c = stokes_drift_profile_windsea_swell(su_surface, sv_surface,
    swell_mean_direction_to=swell_dir_to, swell_mean_period=swell_mean_period,
    swell_height=swell_height,
    wind_sea_mean_direction_to=wind_sea_mean_dir_to, wind_sea_mean_period=wind_sea_mean_period,
    wind_sea_height=wind_sea_height, z=z)

profiles = [
    {'u': su_m, 'v': sv_m, 'z': z, 'kwargs': {'label': 'Stokes, Monochromatic'}},
    {'u': su_e, 'v': sv_e, 'z': z, 'kwargs': {'label': 'Stokes, Exponential'}},
    {'u': su_p, 'v': sv_p, 'z': z, 'kwargs': {'label': 'Stokes, Phillips'}},
    {'u': su_c, 'v': sv_c, 'z': z, 'kwargs': {'label': 'Stokes, Windsea + Swell'}}
    ]

plot_stokes_profile(profiles, view=['vertical', 'birdseye', 'u', 'v'])
