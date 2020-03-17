#!/usr/bin/env python
"""
Openoil 3d
==================================
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil3D import OpenOil3D

o = OpenOil3D(loglevel=20, weathering_model='noaa')

print(o.oiltypes)  # Print available oil types

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([reader_norkyst, reader_arome])

#%%
# Seeding some particles
time = reader_arome.start_time
oiltype = 'GULLFAKS, EXXON'
oiltype = 'ARABIAN MEDIUM, API'
oiltype = 'ALGERIAN CONDENSATE'

#%%
# Seed oil elements at defined position and time
o.seed_elements(lon=4.9, lat=60.1, radius=3000, number=2000,
                time=time, z=0, oiltype=oiltype)

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
#o.plot(filename='openoil3d_drift')
o.plot_vertical_distribution()
o.plot_property('water_fraction', mean=True)
o.plot_property('z')
#o.plot_property('mass_evaporated')
#o.plot_property('water_fraction')
#o.plot_property('interfacial_area')
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_openoil3d_0.gif

