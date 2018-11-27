#!/usr/bin/env python

# Illustrating the difference between Euler and Runge-Kutta propagation
# schemes, using an idealised (analytical) eddy current field.

from datetime import datetime, timedelta
import matplotlib.pyplot as plt

from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift


#########################################
# Setting some parameters
duration = timedelta(seconds=15)  # T
time_step=timedelta(seconds=.5)
time_step_output=timedelta(seconds=.5)
delta=.01  # spatial resolution
steps = int(duration.total_seconds()/
            time_step_output.total_seconds() + 1)
#########################################

o = OceanDrift(loglevel=20)

o.fallback_values['land_binary_mask'] = 0
# Note that Runge-Kutta here makes a difference to Euler scheme
o.set_config('drift:scheme', 'runge-kutta4')

double_gyre = reader_double_gyre.Reader(epsilon=.25, omega=0.628, A=0.1)
print double_gyre

o.add_reader(double_gyre)

# Calculate Lyapunov exponents
times = [double_gyre.initial_time +
         n*time_step_output for n in range(steps)]
lcs = o.calculate_ftle(time=times, time_step=time_step,
                       duration=duration, delta=delta, RLCS=False)

# Make run with particles for the same period
o.reset()
x = [.9]
y = [.5]
lon, lat = double_gyre.xy2lonlat(x, y)

o.seed_elements(lon, lat, radius=.1, number=1000,
                time=double_gyre.initial_time)
o.run(duration=duration, time_step=time_step,
      time_step_output=time_step_output)
o.animation(buffer=0, lcs=lcs, filename='ALCS_with_particles.mp4')

