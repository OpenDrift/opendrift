#!/usr/bin/env python
"""
Entrainment rate
==================================
"""

from datetime import datetime, timedelta
from opendrift.models.openoil3D import OpenOil3D
import matplotlib.pyplot as plt

####################################################
# Tkcalich & Chan (2002) entrainment rate
####################################################
o = OpenOil3D(loglevel=0, weathering_model='noaa')
o.fallback_values['land_binary_mask'] = 0
o.fallback_values['x_sea_water_velocity'] = -.2
o.fallback_values['y_sea_water_velocity'] = 0
o.fallback_values['x_wind'] = 12
o.fallback_values['y_wind'] = 0
o.fallback_values['sea_surface_wave_stokes_drift_x_velocity'] = .3
o.fallback_values['sea_surface_wave_stokes_drift_y_velocity'] = 0
o.set_config('wave_entrainment:entrainment_rate', 'Tkalich & Chan (2002)')
o.set_config('wave_entrainment:droplet_size_distribution', 'Johansen et al. (2015)')
o.set_config('processes:evaporation', False)
o.set_config('processes:dispersion', False)
o.set_config('turbulentmixing:droplet_diameter_min_wavebreaking', 1e-6)
o.set_config('turbulentmixing:droplet_diameter_max_wavebreaking', 1e-3)
o.seed_elements(lon=4, lat=60, time=datetime.now(), number=1000,
                radius=100, z=0, oiltype='VILJE')
o.run(duration=timedelta(hours=24), time_step=900)

######################################################
# Li et al. (2017) entrainment rate
######################################################
o2 = OpenOil3D(loglevel=0, weathering_model='noaa')
o2.fallback_values['land_binary_mask'] = 0
o2.fallback_values['x_sea_water_velocity'] = -.2
o2.fallback_values['y_sea_water_velocity'] = 0
o2.fallback_values['x_wind'] = 12
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
                radius=100, z=0, oiltype='VILJE')
o2.run(duration=timedelta(hours=24), time_step=900)

###########################
# Plotting and comparing
###########################
o.plot_vertical_distribution()
o2.plot_vertical_distribution()
o.plot_oil_budget()
o2.plot_oil_budget()
legend = ['Tkalich & Chan (2002)', 'Li et al. (2017)']
o.animation_profile(compare=o2, legend=legend, filename='entrainment_rate.gif')

#%%
# .. image:: /gallery/animations/entrainment_rate.gif

o.animation(compare=o2, legend=legend, filename='entrainment_rate_compare.gif')


#%%
# .. image:: /gallery/animations/entrainment_rate_compare.gif



