#!/usr/bin/env python

from datetime import datetime, timedelta
from opendrift.models.openoil3D import OpenOil3D

import matplotlib.pyplot as plt
import numpy as np

######################################################
# Li et al. (2017) entrainment rate (light vs. heavy oil)
######################################################
o2 = OpenOil3D(loglevel=0, weathering_model='noaa')
o2.fallback_values['land_binary_mask'] = 0
o2.fallback_values['x_sea_water_velocity'] = -.2
o2.fallback_values['y_sea_water_velocity'] = 0
o2.fallback_values['x_wind'] = 10
o2.fallback_values['y_wind'] = 0
o2.fallback_values['sea_surface_wave_stokes_drift_x_velocity'] = .3
o2.fallback_values['sea_surface_wave_stokes_drift_y_velocity'] = 0
o2.set_config('wave_entrainment:entrainment_rate', 'Li et al. (2017)')
o2.set_config('wave_entrainment:droplet_size_distribution', 'Johansen et al. (2015)')
o2.set_config('processes:evaporation', False)
o2.set_config('processes:dispersion', False)
o2.set_config('turbulentmixing:droplet_diameter_min_wavebreaking', 1e-6)
o2.set_config('turbulentmixing:droplet_diameter_max_wavebreaking', 1e-3)
o2.seed_elements(lon=4, lat=60, time=datetime.now(), number=1000,
                radius=100, z=0, oiltype='TIA JUANA HEAVY, OIL & GAS')
o2.run(duration=timedelta(hours=12), time_step=900)

o3 = OpenOil3D(loglevel=0, weathering_model='noaa')
o3.fallback_values['land_binary_mask'] = 0
o3.fallback_values['x_sea_water_velocity'] = -.2
o3.fallback_values['y_sea_water_velocity'] = 0
o3.fallback_values['x_wind'] = 10
o3.fallback_values['y_wind'] = 0
o3.fallback_values['sea_surface_wave_stokes_drift_x_velocity'] = .3
o3.fallback_values['sea_surface_wave_stokes_drift_y_velocity'] = 0
o3.set_config('wave_entrainment:entrainment_rate', 'Li et al. (2017)')
o3.set_config('wave_entrainment:droplet_size_distribution', 'Johansen et al. (2015)')
o3.set_config('processes:evaporation', False)
o3.set_config('processes:dispersion', False)
o3.set_config('turbulentmixing:droplet_diameter_min_wavebreaking', 1e-6)
o3.set_config('turbulentmixing:droplet_diameter_max_wavebreaking', 1e-3)
o3.seed_elements(lon=4, lat=60, time=datetime.now(), number=1000,
                 radius=100, z=0, oiltype='TIA JUANA LIGHT, OIL & GAS') #'EKOFISK BLEND, STATOIL' similar ent.
o3.run(duration=timedelta(hours=12), time_step=900)

###########################
# Plotting and comparing
###########################
print('#######################')
print('Entrainment rate (heavy)', np.mean(o2.oil_wave_entrainment_rate()))
print('Entrainment rate (light)', np.mean(o3.oil_wave_entrainment_rate()))
print('Viscosity (heavy)', np.mean(o2.elements.viscosity))
print('Viscosity (light)', np.mean(o3.elements.viscosity))
print('Density (heavy)', np.mean(o2.elements.density))
print('Density (light)', np.mean(o3.elements.density))
print('#######################')

o2.plot_oil_budget()
o3.plot_oil_budget()
legend = ['TIA JUANA HEAVY', 'TIA JUANA LIGHT']
o2.animation_profile(compare=o3, legend=legend)
o2.animation(compare=o3, legend=legend)



