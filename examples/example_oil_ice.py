#!/usr/bin/env python

from datetime import datetime, timedelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=0)
o.max_speed = 1

# Nordc4
reader_arctic = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '2Feb2016_Nordic_sigma_3d/Arctic20_1to5Feb_2016.nc')
reader_arctic.interpolation = 'linearND'

#o.add_reader([reader_basemap, reader_arctic])
#o.add_reader([reader_arctic, reader_basemap])
o.add_reader([reader_arctic])

# Seeding some particles
lon = 26.4; lat = 77.3;  # Spitzbergen

#time = datetime(2015, 9, 22, 6, 0, 0)
#time = [reader_nordic4.start_time,
#        reader_nordic4.start_time + timedelta(hours=30)]
time = reader_arctic.start_time

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=7000, number=3000, time=time)
o.fallback_values['y_wind'] = 4  # Adding some northwards wind

print(o)

# Adjusting some configuration
o.set_config('general:basemap_resolution',  'i')
o.set_config('processes:dispersion',  False)
o.set_config('processes:evaporation',  False)
o.set_config('processes:emulsification',  False)
o.set_config('drift:current_uncertainty',  .5)
o.set_config('drift:wind_uncertainty',  5)

# Running model (until end of driver data)
o.run(duration=timedelta(days=4), time_step=3600)
print(o.readers['basemap_landmask'])

# Print and plot results
print(o)
o.animation()
o.plot(background='sea_ice_area_fraction')
