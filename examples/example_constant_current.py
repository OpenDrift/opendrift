#!/usr/bin/env python

from datetime import datetime

from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Adding no input models, but instead constant northwards current of 1 m/s
o.fallback_values['x_sea_water_velocity'] = 0
o.fallback_values['y_sea_water_velocity'] = 1
o.fallback_values['land_binary_mask'] = 0

# Seeding some particles
lon = 4.3; lat = 60.0; # Outside Bergen
time = datetime(2015, 9, 22, 6, 0, 0)
# Seed elements at defined position and time
o.seed_elements(lon, lat, radius=5000, number=100, time=time)

print(o)

# Running model for 50 hours
o.run(steps=50)

# Print and plot results
print(o)
o.plot()
