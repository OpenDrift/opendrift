#!/usr/bin/env python
"""
Relative
==================================
"""

from datetime import timedelta

from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

# Landmask
reader_landmask = reader_global_landmask.Reader(
                    llcrnrlon=3, llcrnrlat=59.5,
                    urcrnrlon=7, urcrnrlat=62)

o.add_reader([reader_landmask, reader_norkyst, reader_arome])

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
o.set_config('drift:wind_uncertainty',  1)
o.set_config('drift:relative_wind',  False)

# Running model
o.run(steps=66*2, time_step=1800)

# Second run, for comparison
o2 = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information
o2.add_reader([reader_landmask, reader_norkyst, reader_arome])
o2.seed_elements(lon, lat, radius=50, number=5000, time=time)
o2.set_config('processes:dispersion',  False)
o2.set_config('processes:evaporation',  False)
o2.set_config('processes:emulsification',  False)
o2.set_config('drift:current_uncertainty',  .1)
o2.set_config('drift:wind_uncertainty',  1)
o2.set_config('drift:relative_wind',  True)
o2.run(steps=66*2, time_step=1800)


# Animate and compare the two runs
o.animation(compare=o2, legend=['Absolute wind', 'Relative wind'], filename='relative.gif')

#%%
# .. image:: /gallery/animations/relative.gif

