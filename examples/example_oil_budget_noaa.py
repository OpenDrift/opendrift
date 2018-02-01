#!/usr/bin/env python

from datetime import datetime
from opendrift.models.openoil3D import OpenOil3D

# NOAA OilLibrary must be installed to run this example
o = OpenOil3D(loglevel=0, weathering_model='noaa')

# Using constand wind and current
o.fallback_values['x_wind'] = 7
o.fallback_values['x_sea_water_velocity'] = .7
o.fallback_values['y_sea_water_velocity'] = .3

# Seeding some particles
lon = 4.7; lat = 60.2; # Outside Bergen
time = datetime.now()
o.seed_elements(lon, lat, radius=3000, number=500, time=time, z=0,
                #oiltype='GULLFAKS, EXXON')
                oiltype='ARABIAN MEDIUM, API')
                #oiltype='ALGERIAN CONDENSATE')

# Adjusting some configuration
o.set_config('processes:evaporation', True)
o.set_config('processes:emulsification', True)
o.set_config('processes:turbulentmixing', True)
o.set_config('turbulentmixing:timestep', 2)

# Running model (until end of driver data)
o.run(steps=4*20, time_step=900,
      outfile='oil_budget.nc')

# Print and plot results
print o
o.plot_oil_budget()
o.plot()
o.animation()
o.plot_property('fraction_evaporated')
o.plot_property('density')
o.plot_property('water_fraction')
o.plot_property('viscosity')
o.plot_property('interfacial_area')
o.plot_property('z')
