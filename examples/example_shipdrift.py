#!/usr/bin/env python

from datetime import datetime, timedelta

from opendrift.models.shipdrift import ShipDrift

o = ShipDrift(loglevel=0)
o.set_config('general:basemap_resolution', 'h')

o.add_readers_from_list([
    'http://thredds.met.no/thredds/dodsC/sea/nordic4km/zdepths1h/aggregate_be',
    'http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc',
    'http://thredds.met.no/thredds/dodsC/sea/mywavewam4/mywavewam4_be'
    ])

# Seeding some particles
lon = 5.0; lat = 63.0; # Outside Bergen

#time = reader_arome.start_time
#time = datetime(2016, 12, 20, 17, 0)
time = datetime.now()

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=1000, number=1000, time=time,
                length=80.0, beam=10.0, height=9.0, draft=4.0)

print(o.elements_scheduled)
print(o)

#o.set_projection('+proj=merc')

# Running model
o.run(steps=24, stop_on_error=True)

# Print and plot results
print(o)
o.plot(linecolor='orientation')
#o.animation()
