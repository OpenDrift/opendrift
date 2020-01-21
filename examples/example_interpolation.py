#!/usr/bin/env python
"""
Interpolation
==================================

Comparing simulation with two different interpolation methods:
'ndimage' and 'linearND'. The third option is 'linearNDFast' (default)
"""

from datetime import timedelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([reader_norkyst, reader_arome])

# Seeding some particles
lon = 4.2; lat = 60.0; # Outside Bergen
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=50, number=5000, time=time)


# Adjusting some configuration
o.set_config('processes:dispersion',  False)
o.set_config('processes:evaporation',  False)
o.set_config('processes:emulsification',  False)
o.set_config('drift:current_uncertainty',  .1)
o.set_config('drift:wind_uncertainty',  2)
o.set_config('drift:relative_wind',  False)

for r in o.readers:
    o.readers[r].interpolation = 'ndimage'

# Running model
o.run(steps=66*2, time_step=1800)

# Second run, for comparison
o2 = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information
o2.add_reader([reader_norkyst, reader_arome])
o2.seed_elements(lon, lat, radius=50, number=5000, time=time)
o2.set_config('processes:dispersion',  False)
o2.set_config('processes:evaporation',  False)
o2.set_config('processes:emulsification',  False)
o2.set_config('drift:current_uncertainty',  .1)
o2.set_config('drift:wind_uncertainty',  2)
o2.set_config('drift:relative_wind',  False)

for r in o.readers:
    o.readers[r].interpolation = 'linearND'

o2.run(steps=66*2, time_step=1800)


# Animate and compare the two runs
o.animation(compare=o2, legend=['ndimage', 'linearND'], filename='interpolation.gif')


#%%
# .. image:: /gallery/animations/interpolation.gif
