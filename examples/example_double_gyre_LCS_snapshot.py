#!/usr/bin/env python
"""
Double gyre - LCS with particles
============================================

Drift of particles in an idealised (analytical) eddy current field,
plotted on top of the LCS. This takes some minutes to calculate.
"""

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np

from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift

#%%
# Setting some parameters
plot_time =  timedelta(seconds=12)
duration = timedelta(seconds=12)  # T
time_step=timedelta(seconds=.5)
time_step_output=timedelta(seconds=.5)
delta=.01  # spatial resolution
#steps = int(duration.total_seconds()/ time_step_output.total_seconds() + 1)

o = OceanDrift(loglevel=20)

#%%
# Note that Runge-Kutta here makes a difference to Euler scheme
o.set_config('drift:scheme', 'runge-kutta4')
o.disable_vertical_motion()
o.fallback_values['land_binary_mask'] = 0

double_gyre = reader_double_gyre.Reader(epsilon=.25, omega=0.628, A=0.1)
print(double_gyre)
o.add_reader(double_gyre)


#%%
# Calculate Lyapunov exponents
#times = [double_gyre.initial_time +  n*time_step_output for n in range(steps)]
ftle = o.calculate_ftle(time=double_gyre.initial_time + plot_time, time_step=time_step,
                       duration=duration, delta=delta, RLCS=False)

lcs = o.calculate_lcs(time=double_gyre.initial_time + plot_time, time_step=-time_step,
                       duration=duration, delta=delta)
#%%
# Make run with particles for the same period
o.reset()
x = np.linspace(0.1,1.9,100)
y = np.linspace(0.1,0.9,100)
X,Y = np.meshgrid(x.ravel(),y.ravel())
lon, lat = double_gyre.xy2lonlat(X,Y)

o.seed_elements(lon, lat,
                time=double_gyre.initial_time)
o.disable_vertical_motion()
o.run(duration=plot_time, time_step=time_step,
      time_step_output=time_step_output)
lonmin, latmin = double_gyre.xy2lonlat(0.,0.)
lonmax, latmax = double_gyre.xy2lonlat(2.,1.)
o.plot(lcs=ftle, show_initial=False, linewidth = 0, corners =[lonmin, lonmax, latmin, latmax], cmap='cividis')

fig = plt.figure()
fig.add_subplot(211)
plt.pcolormesh(np.log(np.sqrt(lcs['eigval'][0,:,:,0])), cmap='cividis'),plt.colorbar()
plt.quiver(lcs['eigvec'][0,:,:,0,0], lcs['eigvec'][0,:,:,0,1])
fig.add_subplot(212)
plt.pcolormesh(np.log(np.sqrt(lcs['eigval'][0,:,:,1])), cmap='cividis'),plt.colorbar()
plt.quiver(lcs['eigvec'][0,:,:,0,0], lcs['eigvec'][0,:,:,0,1])

#o.animation(buffer=0, lcs=ftle, hide_landmask=True)

#%%
# .. image:: /gallery/animations/example_double_gyre_LCS_particles_0.gif
