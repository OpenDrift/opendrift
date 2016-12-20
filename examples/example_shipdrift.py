#!/usr/bin/env python

from datetime import datetime, timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.models.shipdrift import ShipDrift

o = ShipDrift(loglevel=0, basemap_resolution='i')

o.add_readers_from_list([
    #'http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be',
    'http://thredds.met.no/thredds/dodsC/fou-hi/nordic4km-1h/Nordic-4km_SURF_1h_avg_00.nc',
    #'http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc',
    '/lustre/storeB/project/metproduction/products/ecmwf/nc/ec_atmo_0_1deg_20161220T120000Z_1h.nc',
    'http://thredds.met.no/thredds/dodsC/sea/mywavewam4/mywavewam4_be'
    ])

# Seeding some particles
lon = 5.0; lat = 63.0; # Outside Bergen

#time = reader_arome.start_time
time = datetime(2016, 12, 20, 17, 0)
#time = datetime.now()

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=1000, number=5, time=time)

print o.elements_scheduled
print o

#o.set_projection('+proj=merc')

# Running model
o.run(steps=24, stop_on_error=True)

# Print and plot results
print o
o.plot(linecolor='orientation')
#o.animation()
