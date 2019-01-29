#!/usr/bin/env python

from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from opendrift.models.openoil3D import OpenOil3D


number = 10000
steps = 20
mass_oil = 2000  # mass oil per particle
oiltype = '*GENERIC DIESEL'

####################################################
# First run, where surface oil thickness is updated 
####################################################
o1 = OpenOil3D(loglevel=30, weathering_model='noaa')
# Northwards wind, eastwards current
o1.fallback_values['land_binary_mask'] = 0
o1.fallback_values['x_wind'] = 0
o1.fallback_values['y_wind'] = 7
o1.fallback_values['sea_surface_wave_stokes_drift_x_velocity'] = 0
o1.fallback_values['sea_surface_wave_stokes_drift_y_velocity'] = .3
o1.fallback_values['x_sea_water_velocity'] = .1
o1.fallback_values['y_sea_water_velocity'] = 0 
# Using Johansen droplet spectrum, which depends on oil film thickness
o1.set_config('wave_entrainment:droplet_size_distribution',
             'Johansen et al. (2015)')
o1.set_config('drift:wind_uncertainty', 2)
o1.set_config('drift:current_uncertainty', .1)
o1.set_config('processes:update_oilfilm_thickness', True)

o1.seed_elements(lon=4.5, lat=60, number=number,
                mass_oil=mass_oil, radius=1000,
                oiltype=oiltype,
                time=datetime.now())
o1.run(steps=steps)

# Animation shows how oil thickness evolves,
# and decreases due to evaporation and spreading
unitfactor=1e6  # show film thickness in micrometers
o1.animation(color='oil_film_thickness',
             vmin=1e-7*unitfactor, vmax=1e-4*unitfactor,
             unitfactor=unitfactor, surface_only=True)

###################################################################
# Second run, identical but without updating surface oil thickness
###################################################################
print('Second run...')
o2 = OpenOil3D(loglevel=30, weathering_model='noaa')
o2.fallback_values = o1.fallback_values

o2.set_config('wave_entrainment:droplet_size_distribution',
             'Johansen et al. (2015)')
o2.set_config('drift:wind_uncertainty', 2)
o2.set_config('drift:current_uncertainty', .1)
o2.set_config('processes:update_oilfilm_thickness', False)

o2.seed_elements(lon=4.5, lat=60, number=number,
                mass_oil=mass_oil, radius=1000,
                oiltype=oiltype,
                time=datetime.now())
o2.run(steps=steps)

####################################################
# Comparison plots
####################################################
# We see that with the updated film thickness,
# the droplets are getting gradually smaller
r1 = o1.get_property('diameter')[0]
r2 = o2.get_property('diameter')[0]
plt.plot(np.median(r1*1e6, 1))
plt.plot(np.median(r2*1e6, 1))
plt.legend(['With updated film thickness', 'With constant film thickness'])
plt.xlabel('Time step')
plt.ylabel('Median droplet diameter  [micrometer]')
plt.show()

# We see that oil film thickness has virtually no impact on horizontal drift
o1.animation(compare=o2, legend=['Updated film thickness',
                                 'Constant/default film thickness'])
