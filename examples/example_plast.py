#!/usr/bin/env python

from datetime import datetime, timedelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.plastdrift import PlastDrift

o = PlastDrift(loglevel=0)

# Arome atmospheric model
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
# Norkyst ocean model
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([reader_norkyst, reader_arome])
o.set_config('general:basemap_resolution', 'h')

# Seeding some particles
lon = 4.6; lat = 60.0; # Outside Bergen
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]
o.seed_elements(lon, lat, radius=50, number=3000, time=time)
                
o.run(end_time=reader_norkyst.end_time, time_step=1800,
      time_step_output=3600)

# Print and plot results
print(o)
o.animation(color='z')
o.plot_property('z')
