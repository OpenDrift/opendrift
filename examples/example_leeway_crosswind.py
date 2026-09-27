#!/usr/bin/env python
"""
Leeway crosswind direction and jibing
======================================

Idealised test with constant 10 m/s wind blowing towards north, and no current.
Leeway coefficients are given for drift to the right (CWR) and left (CWL)
of downwind, with positive crosswind leeway to the right of downwind.

PIW-4 (person in survival suit) is asymmetric:

- right: slope 1.36 %, offset -3.30 cm/s -> about 10.3 cm/s towards east
- left: slope -0.13 %, offset -2.65 cm/s -> about 3.95 cm/s towards west

After 24 hours without jibing, right-drifting elements should thus be about 8.9 km
east and left-drifting elements about 3.4 km west of the downwind axis.
"""

from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from opendrift.models.leeway import Leeway


def run(object_type, jibe_probability):
    o = Leeway(loglevel=50)
    o.set_config('environment:constant', {'x_sea_water_velocity': 0, 'y_sea_water_velocity': 0,
                    'x_wind': 0, 'y_wind': 10, 'land_binary_mask': 0})
    o.seed_elements(lon=4, lat=60, time=datetime(2020, 1, 1), number=2000,
                    object_type=object_type, jibe_probability=jibe_probability)
    o.run(duration=timedelta(hours=24), time_step=900)
    x = (o.elements.lon - 4) * 111.2 * np.cos(np.radians(60))
    y = (o.elements.lat - 60) * 111.2
    return o, x, y


#%%
# PIW-4 without jibing: elements seeded as right-drifting (orientation 0)
# should end up east of the downwind axis, and left-drifting west.
o, x, y = run(object_type=4, jibe_probability=0)
ori = o.elements.orientation
for orientation, name, expected in ((0, 'right', 8.9), (1, 'left', -3.4)):
    print(f'{name:5s}: mean crosswind {x[ori == orientation].mean():5.1f} km '
          f'(expected {expected:5.1f} km)')

fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
for ax, jibe_probability in zip(axes, (0, 0.04)):
    o, x, y = run(object_type=4, jibe_probability=jibe_probability)
    ori = o.elements.orientation
    ax.scatter(x[ori == 0], y[ori == 0], s=2, c='tab:red', label='right of downwind')
    ax.scatter(x[ori == 1], y[ori == 1], s=2, c='tab:blue', label='left of downwind')
    ax.plot(0, 0, 'k*', markersize=12)
    ax.annotate('', xy=(0, 5), xytext=(0, 0), arrowprops=dict(arrowstyle='->'))
    ax.set_title(f'PIW-4, 24 h, jibe probability {jibe_probability}/h')
    ax.set_xlabel('Crosswind (east) [km]')
    ax.axvline(0, color='gray', lw=.5)
    ax.set_aspect('equal')
axes[0].set_ylabel('Downwind (north) [km]')
axes[0].legend(loc='lower left', markerscale=5)
plt.show()

#%%
# Spread of the generic PIW-1 category compared to more specific PIW classes.
# PIW-1 is a mean over all persons-in-water, with large standard deviations
# of the leeway coefficients, and a much larger search area results.
fig, ax = plt.subplots(figsize=(7, 7))
for object_type, color in ((1, 'tab:gray'), (3, 'tab:green'), (6, 'tab:orange')):
    o, x, y = run(object_type=object_type, jibe_probability=0.04)
    key = o.leewayprop[object_type]['OBJKEY']
    print(f'{key}: downwind std {y.std():4.1f} km, crosswind std {x.std():4.1f} km')
    ax.scatter(x, y, s=2, c=color, label=key)
ax.plot(0, 0, 'k*', markersize=12)
ax.set_xlabel('Crosswind (east) [km]')
ax.set_ylabel('Downwind (north) [km]')
ax.set_aspect('equal')
ax.legend(markerscale=5)
ax.set_title('10 m/s wind towards north, 24 hours')
plt.show()
