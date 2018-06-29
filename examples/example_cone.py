#!/usr/bin/env python

from datetime import timedelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([reader_norkyst, reader_arome])

# Seeding some particles
lon = [3.6, 5.1]; lat = [61., 59.6]; # Outside Bergen

# Seed elements along cone, e.g. ship track with
# increasing uncertainty in position
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]
#time = reader_arome.start_time

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=[1000, 10000], number=5000,
                time=time, cone=True)

# Adjusting some configuration
o.set_config('processes:dispersion', True)
o.set_config('processes:evaporation', True)
o.set_config('processes:emulsification', True)
o.set_config('drift:current_uncertainty', .1)
o.set_config('drift:wind_uncertainty', 1)

# Running model (until end of driver data)
o.run(steps=66*2, time_step=1800)

# Print and plot results
print(o)
o.plot()
o.animation()
