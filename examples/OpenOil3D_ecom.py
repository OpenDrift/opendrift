#!/usr/bin/env python
"""
Openoil 3d_ECOM
==================================
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_ECOM_building
from opendrift.models.openoil3D import OpenOil3D

o = OpenOil3D(loglevel=0, weathering_model='noaa')

#print(o.oiltypes)  # Print available oil types


reader_pcse = reader_ECOM_building.Reader('/home/arian/IC/arquivos_cdf/z_six_days.nc')

o.add_reader([reader_pcse])

#%%
# Seeding some particles
oiltype = 'GULLFAKS, EXXON'
oiltype = 'ARABIAN MEDIUM, API'
oiltype = 'ALGERIAN CONDENSATE'

#%%
lon = -45.1; lat = -24;

time = reader_pcse.start_time
# Seed elements at defined position and time
import numpy as np

o.seed_elements(lon, lat,  radius=200, number=200,
                time=time,z= -45, oiltype=oiltype)

#%%
# Adjusting some configuration
o.set_config('processes:evaporation',  True)
o.set_config('processes:emulsification',  True)
o.set_config('processes:turbulentmixing',  True)
o.set_config('turbulentmixing:timestep',  5)

#%%
# Running model (until end of driver data)
o.run(steps=4*40, time_step=900,
      time_step_output=3600)

#%%
# Print and plot results
print(o)
o.plot(fast=True)
o.plot_oil_budget()
o.plot(filename='openoil3d_drift')
o.plot_vertical_distribution()
o.plot_property('water_fraction', mean=True)
o.plot_property('z')
o.plot_property('mass_evaporated')
o.plot_property('water_fraction')
o.plot_property('interfacial_area')
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_openoil3d_0.gif

