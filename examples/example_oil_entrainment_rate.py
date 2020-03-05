#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from opendrift.models.physics_methods import oil_wave_entrainment_rate_li2017, oil_wave_entrainment_rate_tkalich2002

# Viscosities from 0 to 20 Stokes
vis = np.linspace(0, 20, 100)

colors = ['b', 'g', 'r']
lstyles = ['-', '-.']

for wind, color in zip([6, 10, 15], colors):

    # Li (2017)
    for ift, ls in zip([.003, .03], lstyles):
        r = oil_wave_entrainment_rate_li2017(
                dynamical_viscosity=vis, oil_density=950, interfacial_tension=ift,
                wind_speed=wind)
        p1h = 1-np.exp(-3600*r)
        plt.plot(vis, p1h, color+ls, label='Li(2017), %s m/s wind, IFT: %s' % (wind, ift))

    # Tkalich (2002)
    r = oil_wave_entrainment_rate_tkalich2002(wind_speed=wind)
    p1h = 1-np.exp(-3600*r)*np.ones(vis.shape)
    plt.plot(vis, p1h, color + '.', label='Tkalich, wind = %s m/s' % wind)

plt.ylabel('Fraction entrained in 1 hour')
plt.xlabel('Dynamical viscosity [Stokes]')
plt.legend()
plt.show()
