#!/usr/bin/env python

from datetime import datetime

from opendrift.models.openoil import OpenOil

try:
    import gdal
except:
    raise ValueError('Please install GDAL (www.gdal.org) with Python'
                     ' support to run this example')

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information
o.max_speed = .5  # To minimise plotting boundaries

#####################################################
# Seed oil particles within contours from shapefile
#####################################################
o.seed_from_shapefile(o.test_data_folder() +
                      'shapefile_spawning_areas/Torsk.shp',
                      number=2000, layername=None,
                      featurenum=[2, 4], time=datetime.now())

o.fallback_values['x_wind'] = -4  # Constant wind drift
o.fallback_values['y_wind'] = 8
o.set_config('drift:wind_uncertainty', 4) # Adding some diffusion

# Running model
o.run(steps=50, time_step=3600)

# Print and plot results
print(o)
o.plot()
o.animation()
