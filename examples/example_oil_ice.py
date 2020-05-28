#!/usr/bin/env python
"""
Oil in ice
==================================
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20)

# Nordc4
reader_arctic = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '2Feb2016_Nordic_sigma_3d/Arctic20_1to5Feb_2016.nc')
o.add_reader([reader_arctic])

#time = datetime(2015, 9, 22, 6, 0, 0)
#time = [reader_nordic4.start_time,
#        reader_nordic4.start_time + timedelta(hours=30)]
time = reader_arctic.start_time

#%%
# Seed oil elements at defined position and time
#o.seed_elements(lon=24.4, lat=77.3, radius=7000, number=3000, time=time)
o.seed_elements(lon=27, lat=77.0, radius=5000, number=3000, time=time)
o.fallback_values['y_wind'] = 7  # Adding some northwards wind

#%%
# Adjusting some configuration
o.set_config('processes:dispersion',  False)
o.set_config('processes:evaporation',  False)
o.set_config('processes:emulsification',  False)
o.set_config('drift:current_uncertainty',  .5)
o.set_config('drift:wind_uncertainty',  3)

#%%
# Running model (until end of driver data)
o.run(duration=timedelta(days=4), time_step=3600, time_step_output=3600*3)

#%%
# Print and plot results
print(o)
o.animation(background='sea_ice_area_fraction', fast=True, buffer=2)

#%%
# .. image:: /gallery/animations/example_oil_ice_0.gif

o.plot(background='sea_ice_area_fraction', fast=True)
