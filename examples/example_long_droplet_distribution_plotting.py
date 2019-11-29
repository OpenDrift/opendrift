#!/usr/bin/env python
"""
Droplet distribution (plotting)
==================================

Plotting different droplet size distributions used in Opendrift (see openoil3D.py script)
"""

import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime, timedelta
from opendrift.models.openoil3D import OpenOil3D

####################################################
# Delvigne & Sweeney (1988) droplet spectrum
####################################################
o = OpenOil3D(loglevel=0, weathering_model='noaa')
o.fallback_values['land_binary_mask'] = 0
o.fallback_values['x_sea_water_velocity'] = -.2
o.fallback_values['y_sea_water_velocity'] = 0
o.fallback_values['x_wind'] = 12
o.fallback_values['y_wind'] = 0
o.fallback_values['sea_surface_wave_stokes_drift_x_velocity'] = .3
o.fallback_values['sea_surface_wave_stokes_drift_y_velocity'] = 0
o.set_config('wave_entrainment:droplet_size_distribution', 'Exponential')
o.set_config('processes:evaporation', False)
o.set_config('processes:dispersion', False)
o.set_config('turbulentmixing:droplet_diameter_min_wavebreaking', 1e-6)
o.set_config('turbulentmixing:droplet_diameter_max_wavebreaking', 1e-3)
o.set_config('turbulentmixing:droplet_size_exponent', -2.3)
o.seed_elements(lon=4, lat=60, time=datetime.now(), number=1000, radius=100,
                z=0, oiltype='VILJE')
o.run(duration=timedelta(hours=1), time_step=3600)

######################################################
# Uniform droplet spectrum
######################################################
o2 = OpenOil3D(loglevel=0, weathering_model='noaa')
o2.fallback_values['land_binary_mask'] = 0
o2.fallback_values['x_sea_water_velocity'] = -.2
o2.fallback_values['y_sea_water_velocity'] = 0
o2.fallback_values['x_wind'] = 12
o2.fallback_values['y_wind'] = 0
o2.fallback_values['sea_surface_wave_stokes_drift_x_velocity'] = .3
o2.fallback_values['sea_surface_wave_stokes_drift_y_velocity'] = 0
o2.set_config('wave_entrainment:droplet_size_distribution', 'Exponential')
o2.set_config('processes:evaporation', False)
o2.set_config('processes:dispersion', False)
o2.set_config('turbulentmixing:droplet_diameter_min_wavebreaking', 1e-6)
o2.set_config('turbulentmixing:droplet_diameter_max_wavebreaking', 2e-3)
o2.set_config('turbulentmixing:droplet_size_exponent', 0)
o2.seed_elements(lon=4, lat=60, time=datetime.now(), number=1000, radius=100,
                 z=0, oiltype='VILJE')
o2.run(duration=timedelta(hours=1), time_step=3600)

########################################################
# Johansen et al. (2015) droplet spectrum
########################################################
o3 = OpenOil3D(loglevel=0, weathering_model='noaa')
o3.fallback_values['land_binary_mask'] = 0
o3.fallback_values['x_sea_water_velocity'] = -.2
o3.fallback_values['y_sea_water_velocity'] = 0
o3.fallback_values['x_wind'] = 12
o3.fallback_values['y_wind'] = 0
o3.fallback_values['sea_surface_wave_stokes_drift_x_velocity'] = .3
o3.fallback_values['sea_surface_wave_stokes_drift_y_velocity'] = 0
o3.set_config('wave_entrainment:droplet_size_distribution', 'Johansen et al. (2015)')
o3.set_config('processes:evaporation', False)
o3.set_config('processes:dispersion', False)
o3.set_config('turbulentmixing:droplet_diameter_min_wavebreaking', 1e-6)
o3.set_config('turbulentmixing:droplet_diameter_max_wavebreaking', 2e-3)
o3.seed_elements(lon=4, lat=60, time=datetime.now(), number=1000, radius=100,
                 z=0, oiltype='VILJE', oil_film_thickness=0.001)
o3.run(duration=timedelta(hours=1), time_step=3600)

dmin = 1e-6 # Change dmin, dmax according to choice in config
dmax = 2e-3
droplet_diameters = o.elements.diameter
droplet_diameters2 = o2.elements.diameter
droplet_diameters3 = o3.elements.diameter

################## Plotting ##########################
plt.subplot(2,2,1)
plt.hist(droplet_diameters, 100, range=(dmin,dmax), align='mid')
plt.xlabel('Droplet diameter d [m]', fontsize=8)
plt.ylabel('N(d)', fontsize=8)
plt.title('Delvigne & Sweeney (1988) distribution', fontsize=10)

plt.subplot(2,2,2)
plt.hist(droplet_diameters2, 100, range=(dmin,dmax), align='mid')
plt.xlabel('Droplet diameter d [m]', fontsize=8)
plt.ylabel('N(d)', fontsize=8)
plt.title('Uniform distribution', fontsize=10)

plt.subplot(2,2,3)
plt.hist(droplet_diameters3, 100, range=(dmin,dmax), align='mid')
plt.xlabel('Droplet diameter d [m]', fontsize=8)
plt.ylabel('N(d)', fontsize=8)
plt.title('Johansen et al. (2015) distribution', fontsize=10)
plt.show()


