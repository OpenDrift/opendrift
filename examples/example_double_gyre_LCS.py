#!/usr/bin/env python
"""
Double gyre - Lagrangian Coherent Structures
============================================

Calculating attracting and repelling LCS for an idealised (analytical) eddy current field.
"""

from datetime import datetime, timedelta
import matplotlib.pyplot as plt

from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

o.set_config('environment:fallback:land_binary_mask', 0)
#%%
# Note that Runge-Kutta here makes a difference to Euler scheme
o.set_config('drift:advection_scheme', 'runge-kutta4')

double_gyre = reader_double_gyre.Reader(epsilon=.25, omega=0.628, A=0.1)
print(double_gyre)

o.add_reader(double_gyre)

lcs = o.calculate_ftle(time=double_gyre.initial_time+timedelta(seconds=3),
                       time_step=timedelta(seconds=.5),
                       duration=timedelta(seconds=15),
                       delta=.02)

#%%
# These plots should reproduce Mov 12 on this page:
# https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html
plt.subplot(2,1,1)
plt.imshow(lcs['RLCS'][0,:,:], interpolation='nearest', cmap='jet', origin='lower')
plt.colorbar()
plt.title('Repelling LCS (forwards)')
plt.subplot(2,1,2)
plt.imshow(lcs['ALCS'][0,:,:], interpolation='nearest', cmap='jet', origin='lower')
plt.colorbar()
plt.title('Attracting LCS (backwards)')
plt.show()
