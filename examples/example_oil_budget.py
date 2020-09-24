#!/usr/bin/env python
"""
Oil budget
==================================
"""

from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20, weathering_model='noaa')

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader([reader_norkyst, reader_arome])
#o.fallback_values['x_wind'] = 9
#o.fallback_values['y_wind'] = 0
#o.fallback_values['x_sea_water_velocity'] = .7
#o.fallback_values['y_sea_water_velocity'] = .3
#o.fallback_values['land_binary_mask'] = 0

#%%
# Seed oil elements at defined position and time
oiltype = 'IFO-380LS 2014'
oiltype = 'IFO 300'
oiltype = 'IFO-180NS 2014'
oiltype = '*GENERIC LIGHT CRUDE'
oiltype = '*GENERIC HEAVY CRUDE'
o.seed_elements(lon=4.8, lat=60.0, z=0, radius=3000, number=1000,
                time=reader_arome.start_time, oiltype=oiltype)

#%%
# Adjusting some configuration
o.set_config('processes:dispersion', True)
o.set_config('processes:evaporation', True)
o.set_config('processes:emulsification', True)
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:timestep', 20.) # seconds

#%%
# Running model
o.run(duration=timedelta(hours=24), time_step=1800)

o.plot_oil_budget(show_density_viscosity=True, show_wind_and_current=True)

o.animation(color='viscosity')

#%%
# .. image:: /gallery/animations/example_oil_budget_0.gif
