#!/usr/bin/env python
"""
Seed demonstration
=============
"""

from datetime import datetime, timedelta
import numpy as np
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=50)

# We do not care about landmask or current for this seeding demonstration,
# so we simple specify fallback_values instead of adding any readers
o.fallback_values['land_binary_mask'] = 0
o.fallback_values['x_sea_water_velocity'] = 0
o.fallback_values['y_sea_water_velocity'] = 0
o.set_config('drift:wind_uncertainty', 0)
o.set_config('drift:current_uncertainty', 0)

time=datetime(2016, 1, 20, 12, 30, 0)


print('\n' + '='*70)
print('Seeding a single element at a point:')
print('o.seed_elements(lon=4, lat=60, time=time)')
print('='*70)
o.seed_elements(lon=4, lat=60, time=time)
#o.run(steps=1)
o.plot(buffer=.7)


print('\n' + '='*70)
print('Seeding 100 elements within a radius of 1000 m:')
print('o.seed_elements(lon=4, lat=60, number=100, radius=1000, time=time)')
print('='*70)
o.seed_elements(lon=4, lat=60, number=100, radius=1000, time=time)
#o.run(steps=1)
o.plot(buffer=.7)


print('\n' + '='*70)
print('Seeding 100 elements within a radius of 1000 m\n and specifying a property (here z/depth) as an array:')
print('o.seed_elements(lon=4, lat=60, number=100, radius=1000, time=time, z=z)')
print('='*70)
z = np.linspace(0, -50, 100)  # Linearly increasing depth
o.seed_elements(lon=4, lat=60, number=100, radius=1000, time=time, z=z)
o.fallback_values['y_sea_water_velocity'] = 3  # Adding some current to be able to visualise depth as color of trajectories
o.run(steps=1)
o.plot(linecolor='z', buffer=.7)
o.fallback_values['y_sea_water_velocity'] = 0


print('\n' + '='*70)
print('Seeding 100 elements at user defined locations\n (here along line between two points):')
print('lats = np.linspace(60, 61, 100)\n' \
      'lons = np.linspace(4, 4.8, 100)\n' \
      'o.seed_elements(lon=lons, lat=lats, time=time)')
print('='*70)
lats = np.linspace(60, 61, 100)
lons = np.linspace(4, 4.8, 100)
o.seed_elements(lon=lons, lat=lats, time=time)
o.run(steps=1)
o.plot(buffer=.2)


print('\n' + '='*70)
print('Seeding 100 elements between two points specifying keyword "cone=True"\n (achieving the same as previous example):')
print('o.seed_elements(lon=[4, 4.8], lat=[60, 61], number=100, cone=True, time=time)')
print('='*70)
o.seed_elements(lon=[4, 4.8], lat=[60, 61], number=100, cone=True, time=time)
o.run(steps=1)
o.plot(buffer=.2)


print('\n' + '='*70)
print('Seeding 1000 elements along cone with radius/uncertainty\n increasing linearly from 0 to 5000 m:')
print('o.seed_elements(lon=[4, 4.8], lat=[60, 61], number=1000, radius=[0, 5000], cone=True, time=time)')
print('='*70)
o.seed_elements(lon=[4, 4.8], lat=[60, 61], number=1000, radius=[0, 5000], cone=True, time=time)
o.run(steps=1)
o.plot(buffer=.2)


print('\n' + '='*70)
print('If specifying time as a two element list (start and end,\n here +5 hours), elements are seeded linearly in time:')
print('o.seed_elements(lon=[4, 4.8], lat=[60, 61], number=1000, radius=[0, 5000], cone=True, time=[time, time+timedelta(hours=5)])')
print('='*70)
o.seed_elements(lon=[4, 4.8], lat=[60, 61], number=1000, radius=[0, 5000], cone=True, time=[time, time+timedelta(hours=5)])
o.run(steps=5*4, time_step=900)
o.animation(filename='seed_demonstration.gif')

#%%
# .. image:: /gallery/animations/seed_demonstration.gif


print('\n' + '='*70)
print('Any model/module may provide specialised seeding-functions, such as \n seeding oil within contours read from a GML file:')
print('o.seed_from_gml(o.test_data_folder() + "radarsat_oil_satellite_observation/RS2_20151116_002619_0127_SCNB_HH_SGF_433012_9730_12182143_Oil.gml", num_elements=2000)')
print('='*70)
o.reset()
o.seed_from_gml(o.test_data_folder() + 'radarsat_oil_satellite_observation/RS2_20151116_002619_0127_SCNB_HH_SGF_433012_9730_12182143_Oil.gml', num_elements=2000)
o.run(steps=1)
o.plot(buffer=.03)
