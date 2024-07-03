#!/usr/bin/env python
"""
Coastline interaction
======================

Example to illustrate stranding options using an artificial
east-west oscillating current field
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_oscillating
from opendrift.models.oceandrift import OceanDrift
duration = timedelta(hours=36)
time_step = timedelta(hours=.5)
number = 500

# Oscillating east-west current to illustrate stranding options
reader_osc = reader_oscillating.Reader('x_sea_water_velocity', amplitude=1,
                                       zero_time=datetime.utcnow())

#%%
# Coastline option "stranding"
# ============================
#
# Particles are moved forward with forcing velocity (e.g. current and wind drift) times the time step.
# If the next position is on land, the particle is deactivated.
# Note that particles may "jump" a distance onto land, depending on the calculation timestep.
o = OceanDrift(loglevel=50)
o.add_reader([reader_osc])
o.set_config('environment:fallback:y_sea_water_velocity', .2)
o.set_config('general:coastline_action', 'stranding')
o.set_config('general:coastline_approximation_precision', None)
o.seed_elements(lon=12.2, lat=67.7, radius=5000, number=number, time=reader_osc.zero_time)
o.run(duration=duration, time_step=time_step)
print(f'Calculation time: {o.timing["total time"]}')
o.animation()

#%%
# .. image:: /gallery/animations/example_coastline_options_0.gif

#%%
# Coastline option "stranding" with higher precision
# ==================================================
#
# By setting config "general:coastline_approximation_precision" to desired accuracy in degrees (default is 0.01 as in this example),
# a more exact coastline crossing is calculated by the deactivated particles.
# Note that with a (too) large compuation time step, particles may still "jump" over islands.
# An alternative to avoid this possibility is to use a smaller timestep for the simulation, though at a larger computational cost.
o = OceanDrift(loglevel=50)
o.add_reader([reader_osc])
o.set_config('environment:fallback:y_sea_water_velocity', .2)
o.set_config('general:coastline_action', 'stranding')
o.set_config('general:coastline_approximation_precision', .001)  # approx 100m
o.seed_elements(lon=12.2, lat=67.7, radius=5000, number=number, time=reader_osc.zero_time)
o.run(duration=duration, time_step=time_step)
print(f'Calculation time: {o.timing["total time"]}')
o.animation()

#%%
# .. image:: /gallery/animations/example_coastline_options_1.gif

#%%
# Coastline option "previous"
# ===========================
#
# Particles hitting land are moved back to previous position in water,
# and may move offshore if currents/winds/forcing change direction.
# Here we see that several particles are jumping over small island.
# Reducing timestep to e.g. 15 minutes will reduce this problem.
o = OceanDrift(loglevel=50)
o.add_reader([reader_osc])
o.set_config('environment:fallback:y_sea_water_velocity', .2)
o.set_config('general:coastline_action', 'previous')
o.seed_elements(lon=12.2, lat=67.7, radius=5000, number=number, time=reader_osc.zero_time)
o.run(duration=duration, time_step=time_step)
print(f'Calculation time: {o.timing["total time"]}')
o.animation()

#%%
# .. image:: /gallery/animations/example_coastline_options_2.gif

#%%
# Coastline option "none"
# =======================
#
# No interaction with land, as used by e.g. WindBlow model (atmospheric drift) 
o = OceanDrift(loglevel=50)
o.add_reader([reader_osc])
o.set_config('environment:fallback:y_sea_water_velocity', .2)
o.set_config('general:coastline_action', 'none')
o.seed_elements(lon=12.2, lat=67.7, radius=5000, number=number, time=reader_osc.zero_time)
o.run(duration=duration, time_step=time_step)
print(f'Calculation time: {o.timing["total time"]}')
o.animation()

#%%
# .. image:: /gallery/animations/example_coastline_options_3.gif
