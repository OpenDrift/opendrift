#!/usr/bin/env python

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift3D import OceanDrift3D

o = OceanDrift3D(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader([reader_norkyst, reader_arome])

# Seeding some particles
lon = 4.8; lat = 60.0; # Outside Bergen

time = None
time = reader_arome.start_time

# Seed elements at defined position and time
import numpy as np
z = -np.random.rand(1000)*50  # Giving elements a random depth
o.seed_elements(lon, lat, z=z, radius=0, number=1000, time=time)

print(o)

# Adjusting some configuration
o.set_config('drift:current_uncertainty', .1)
o.set_config('drift:wind_uncertainty', 2)

# Running model (until end of driver data)
o.run(steps=66*2, time_step=1800, outfile='openoil.nc',
      export_buffer_length=5)  # Writing to netCDF file every 5 time steps

# Print and plot results
print(o)
o.plot(linecolor='z')  # Color lines according to depth
o.animation()
