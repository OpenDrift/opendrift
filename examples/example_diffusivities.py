#!/usr/bin/env python
"""
Vertical diffusivity
===================
"""

#%%
# Plot the depth dependence of vertical diffusivity from the various analytical methods

import numpy as np
import matplotlib.pyplot as plt
from opendrift.models.physics_methods import verticaldiffusivity_Sundby1983, verticaldiffusivity_Large1994, verticaldiffusivity_stepfunction

depth = np.linspace(0, 60, 100)
windspeed = np.arange(0, 20, 5)
colors = ['r', 'g', 'b', 'k']

for w, c in zip(windspeed, colors):
    plt.plot(depth, np.ones(depth.shape)*
                        verticaldiffusivity_Sundby1983(w),
             c + '-', label='Sundby, wind = %s' % w)
    plt.plot(depth, verticaldiffusivity_Large1994(w, depth, 50),
             c + '--', label='Large, wind = %s, MLD=50' % w)
    plt.plot(depth, verticaldiffusivity_Large1994(w, depth, 20),
             c + '-.', label='Large, wind = %s, MLD=20' % w)

plt.plot(depth, verticaldiffusivity_stepfunction(depth),
         '-m', label='Stepfunction')

plt.xlabel('Depth [m]')
plt.ylabel('Vertical diffusivity [m/s2]')
plt.gca().set_xlim([0, depth.max()])
plt.gca().set_ylim([0, None])
plt.legend()
plt.show()

