#!/usr/bin/env python
"""
Oil film thickness
==================================
"""

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
from opendrift.models.openoil import OpenOil


number = 10000
timestep = timedelta(minutes=10)
timestep_output = timedelta(minutes=60)
duration = timedelta(hours=20)
mass_oil = 2000  # mass oil per particle
oil_type = 'GENERIC DIESEL'
#oil_type = 'GENERIC BUNKER C'

#%%
# First run, where surface oil thickness is updated
o1 = OpenOil(loglevel=20, weathering_model='noaa')
#%%
# Northwards wind, eastwards current
o1.set_config('environment:fallback:land_binary_mask', 0)
o1.set_config('environment:fallback:x_wind', 0)
o1.set_config('environment:fallback:y_wind', 7)
o1.set_config('environment:fallback:sea_surface_wave_stokes_drift_x_velocity', 0)
o1.set_config('environment:fallback:sea_surface_wave_stokes_drift_y_velocity', .3)
o1.set_config('environment:fallback:x_sea_water_velocity', .1)
o1.set_config('environment:fallback:y_sea_water_velocity', 0)
#%%
# Using Johansen droplet spectrum, which depends on oil film thickness
o1.set_config('wave_entrainment:droplet_size_distribution',
             'Johansen et al. (2015)')
o1.set_config('drift:wind_uncertainty', 2)
o1.set_config('drift:current_uncertainty', .1)
o1.set_config('processes:dispersion', False)
o1.set_config('processes:update_oilfilm_thickness', True)

o1.seed_elements(lon=4.5, lat=60, number=number,
                mass_oil=mass_oil, radius=1000,
                oil_type=oil_type,
                time=datetime.utcnow())
o1.run(time_step=timestep, time_step_output=timestep_output,
       duration=duration)

#%%
# Animation shows how oil thickness evolves,
# and decreases due to evaporation and spreading
unitfactor=1e6  # show film thickness in micrometers
o1.animation(color='oil_film_thickness', fast=True,
             vmin=1e-7*unitfactor, vmax=1e-4*unitfactor,
             unitfactor=unitfactor, surface_only=True)

#%%
# .. image:: /gallery/animations/example_oil_thickness_0.gif

#%%
# Second run, identical but without updating surface oil thickness
o2 = OpenOil(loglevel=20, weathering_model='noaa')
o2.set_config('environment:fallback:land_binary_mask', 0)
o2.set_config('environment:fallback:x_wind', 0)
o2.set_config('environment:fallback:y_wind', 7)
o2.set_config('environment:fallback:sea_surface_wave_stokes_drift_x_velocity', 0)
o2.set_config('environment:fallback:sea_surface_wave_stokes_drift_y_velocity', .3)
o2.set_config('environment:fallback:x_sea_water_velocity', .1)
o2.set_config('environment:fallback:y_sea_water_velocity', 0)

o2.set_config('wave_entrainment:droplet_size_distribution',
             'Johansen et al. (2015)')
o2.set_config('drift:wind_uncertainty', 2)
o2.set_config('drift:current_uncertainty', .1)
o2.set_config('processes:dispersion', False)
o2.set_config('processes:update_oilfilm_thickness', False)

o2.seed_elements(lon=4.5, lat=60, number=number,
                mass_oil=mass_oil, radius=1000,
                oil_type=oil_type,
                time=datetime.utcnow())
o2.run(time_step=timestep, time_step_output=timestep_output,
       duration=duration)

#%%
# Comparison plots
o1.plot_oil_budget()
o2.plot_oil_budget()

#%%
# Entrainment
b1 = o1.get_oil_budget()
b2 = o2.get_oil_budget()
plt.plot(b1['mass_surface'], '-r', linewidth=2,
            label='Surface, updated thickness')
plt.plot(b1['mass_submerged'], '--r', linewidth=2,
            label='Submerged, updated thickness')
plt.plot(b1['mass_evaporated'], '-.r', linewidth=2,
            label='Evaporated, updated thickness')
plt.plot(b2['mass_surface'], '-b', linewidth=2,
            label='Surface, constant thickness')
plt.plot(b2['mass_submerged'], '--b', linewidth=2,
            label='Submerged, constant thickness')
plt.plot(b2['mass_evaporated'], '-.b', linewidth=2,
            label='Evaporated, constant thickness')
plt.legend()
plt.xlabel('Time step')
plt.show()

#%%
# We see that with the updated film thickness,
# the droplets are getting gradually smaller
plt.plot(1e6*o1.result.diameter.median(dim='trajectory'))
plt.plot(1e6*o2.result.diameter.median(dim='trajectory'))
plt.legend(['With updated film thickness', 'With constant film thickness'])
plt.xlabel('Time step')
plt.ylabel('Median droplet diameter  [micrometer]')
plt.show()

#%%
# We see that oil film thickness has virtually no impact on horizontal drift
o1.animation(compare=o2, fast=True,
             legend=['Updated film thickness',
                     'Constant/default film thickness'])

#%%
# .. image:: /gallery/animations/example_oil_thickness_1.gif

