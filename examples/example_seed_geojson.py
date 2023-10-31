#!/usr/bin/env python
"""
Seeding from GeoJSON string
===========================
"""

from datetime import datetime, timedelta
from opendrift.models.leeway import Leeway
from opendrift.models.openoil import OpenOil

#%%
# Polygon
#--------

o = OpenOil(loglevel=50)
for var in ['x_wind', 'y_wind', 'x_sea_water_velocity', 'y_sea_water_velocity']:
    o.set_config('environment:constant:' + var, 0)
o.seed_from_geojson("""{
      "type": "Feature",
      "geometry": {
        "type": "Polygon",
        "coordinates": [
          [
            [4.0, 60.0], [4.5, 60.0], [4.7, 60.1],
            [4.2, 60.1], [4.0, 60.0]
          ]
        ]
      },
      "properties": {
        "time": "2020-11-06T12:30:00Z",
        "number": 1000,
        "oil_type": "DVALIN 2020",
        "m3_per_hour": 50
      }
    }""")
o.plot(fast=True)


#%%
# Point release at seafloor
#--------------------------

o = OpenOil(loglevel=50)
o.set_config('environment:constant:sea_floor_depth_below_sea_level', 200)
for var in ['x_wind', 'y_wind', 'x_sea_water_velocity', 'y_sea_water_velocity']:
    o.set_config('environment:constant:' + var, 0)
o.seed_from_geojson("""{
      "type": "Feature",
      "geometry": {
        "type": "Point",
        "coordinates": [ 4.0, 60.0 ]
      },
      "properties": {
        "time": ["2020-11-06T12:30:00Z", "2020-11-06T18:30:00Z"],
        "number": 3000,
        "z": "seafloor"
      }
    }""")

o.run(duration=timedelta(hours=6), time_step=300)
o.animation_profile()

#%%
# .. image:: /gallery/animations/example_seed_geojson_0.gif


#%%
# Cone 
#-----
# from (position1, radius1, time1) to (position2, radius2, time2)

o = Leeway(loglevel=50)
for var in ['x_wind', 'y_wind', 'x_sea_water_velocity', 'y_sea_water_velocity']:
    o.set_config('environment:constant:' + var, 0)
o.seed_from_geojson("""{
      "type": "Feature",
      "geometry": {
        "type": "LineString",
        "coordinates": [
            [4.0, 60.0], [4.5, 60.1]
        ]
      },
      "properties": {
        "time": ["2020-11-06T12:30:00Z", "2020-11-06T18:30:00Z"],
        "radius": [0, 2000],
        "number": 3000
      }
    }""")

o.run(duration=timedelta(hours=6))
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_seed_geojson_1.gif
