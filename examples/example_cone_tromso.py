#!/usr/bin/env python

from datetime import timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Using live data from Thredds
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc')
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader([reader_norkyst, reader_arome])

# Seed elements along cone, e.g. ship track with
# increasing uncertainty in position
latstart = 68.988911 
lonstart = 16.040701
#latstart = 69.477754
#lonstart = 16.441702
latend = 69.991446
lonend = 17.760061
lon = [lonstart, lonend]; lat = [latstart, latend]; # Outside Tromso

time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]
o.seed_elements(lon, lat, radius=[100, 500], number=10000,
                time=time, cone=True)

print(o)

# Adjusting some configuration
o.set_config('processes:dispersion', True)
o.set_config('processes:evaporation', False)
o.set_config('processes:emulsification', True)

# Running model (until end of driver data)
o.run(steps=66*2, time_step=1800)

# Print and plot results
print(o)
o.plot()
o.animation()
#o.animation(filename='oilspill_tromsoe.mp4')
