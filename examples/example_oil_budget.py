#!/usr/bin/env python
"""
Oil budget
==================================
"""

from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil3D import OpenOil3D

o = OpenOil3D(loglevel=20)

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
#    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

#o.add_reader([reader_norkyst, reader_arome])
o.fallback_values['x_wind'] = 7
o.fallback_values['y_wind'] = 0
o.fallback_values['x_sea_water_velocity'] = .7
o.fallback_values['y_sea_water_velocity'] = .3
#o.fallback_values['land_binary_mask'] = 0
#o.add_reader([reader_landmask, reader_norkyst])

#%%
# Seed oil elements at defined position and time
o.seed_elements(lon=4.8, lat=60.0, z=0, radius=3000, number=1000,
                time=reader_arome.start_time, oiltype='EKOFISK')

#%%
# Adjusting some configuration
o.set_config('processes:dispersion', True)
o.set_config('processes:evaporation', False)
o.set_config('processes:emulsification', True)
o.set_config('processes:turbulentmixing', True)
o.set_config('turbulentmixing:TSprofiles', False)

# TODO: Sundby scheme does not work here
#o.set_config('turbulentmixing:diffusivitymodel', 'windspeed_Sundby1983')
o.set_config('turbulentmixing:timestep', 2.) # seconds

#%%
# Running model (until end of driver data)
o.run(steps=4*20, time_step=900, export_buffer_length=10,
      outfile='oil_budget.nc')

# Print and plot results
print(o)
o.plot_oil_budget()
o.plot_property('water_fraction')
o.plot_property('water_fraction', mean=True)
o.plot(fast=True)
o.plot_property('mass_oil')
o.plot_property('z')
o.plot_property('mass_evaporated')
o.plot_property('water_fraction')
o.plot_property('interfacial_area')
o.animation()

#%%
# .. image:: /gallery/animations/example_oil_budget_0.gif
