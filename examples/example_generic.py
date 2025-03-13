#!/usr/bin/env python
"""
Generic example
===============
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

# Arome atmospheric model
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
# Norkyst ocean model
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

# Uncomment to use live data from thredds
#reader_arome = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader([reader_norkyst, reader_arome])

#%%
# Seeding some particles
#time = datetime(2015, 9, 22, 6, 0, 0)
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]
#time = reader_arome.start_time

#%%
# Adjusting some configuration
o.set_config('drift:vertical_mixing', False)
o.set_config('processes:dispersion', False)
o.set_config('processes:evaporation', False)
o.set_config('processes:emulsification', True)
o.set_config('drift:current_uncertainty', .1)
o.set_config('drift:wind_uncertainty', 1)

# Seed oil elements at defined position and time
o.seed_elements(lon=4.6, lat=60.0, radius=50, number=3000, time=time,
                wind_drift_factor=.02)

#%%
# Running model
o.run(end_time=reader_norkyst.end_time, time_step=1800,
      time_step_output=3600, outfile='openoil.nc')

#%%
# Print and plot results
o.plot(fast=True)
#o.plot(background=['x_sea_water_velocity', 'y_sea_water_velocity'], buffer=.5)
#o.animation(fast=True)
#o.animation(density=True, show_elements=False, fast=True)

#%%
# Or an animation can be generated with:
o.animation(fast=True)


#%%
# .. image:: /gallery/animations/example_generic_0.gif
