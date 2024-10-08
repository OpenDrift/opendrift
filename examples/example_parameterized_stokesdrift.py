#!/usr/bin/env python
"""
Parameterised Stokesdrift
=========================
"""

#%%
# If surface Stokes drift is not available from a wave model, there are two alternatives:
#    - one can increase the wind_drift_factor by e.g. 1.5%, as the Stokes Drift is typically 1.5% of the wind speed
#    - or the surface Stokes drift can be parameterized from wind speed and fetch distance

#%%
# The latter option is activated in OpenDrift with the config setting
# `o.set_config('drift:use_tabularised_stokes_drift', True)`

#%%
# This activates a paramterisation of Stokes drift with the following method, as implemented by Petter Nygren from SMHI:
# https://opendrift.github.io/_modules/opendrift/models/physics_methods.html#wave_stokes_drift_parameterised
# The code and corresponding plot below shows how the Stokes drift factor (fraction of wind speed) varies with wind speed and fetch (3 different tabulated fetch distances).

import numpy as np
import matplotlib.pyplot as plt
from opendrift.models.physics_methods import wave_stokes_drift_parameterised

for fetch in ['5000', '25000', '50000']:
    wind = ([np.arange(1, 35), np.array([0])])
    sx, sy = wave_stokes_drift_parameterised(wind=wind, fetch=fetch)
    plt.plot(wind[0], sx/wind[0], label=f'Fetch: {fetch[:-3]} km')
plt.xlabel('Wind speed  [m/s]')
plt.ylabel('Stokes drift / wind speed  [ratio]')
plt.legend()
plt.show()
