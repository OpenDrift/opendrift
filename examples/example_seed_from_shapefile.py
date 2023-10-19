#!/usr/bin/env python
"""
Seeding from shapefile
==================================
"""

from datetime import datetime
from opendrift.models.oceandrift import OceanDrift


o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o.set_config('environment:fallback:x_wind', -4)  # Constant wind drift
o.set_config('environment:fallback:y_wind', 8)
o.set_config('drift:wind_uncertainty', 4) # Adding some diffusion
o.set_config('drift:vertical_mixing', False)

#%%
# Seed particles within contours from shapefile
o.seed_from_shapefile(o.test_data_folder() +
                      'shapefile_spawning_areas/Torsk.shp',
                      number=2000, layername=None,
                      featurenum=[2, 4], time=datetime.utcnow())


#%%
# Running model
o.run(steps=50, time_step=3600)

#%%
# Print and plot results
print(o)
o.plot(fast=True)
o.animation(fast=True)

#%%
# .. image:: /gallery/animations/example_seed_from_shapefile_0.gif
