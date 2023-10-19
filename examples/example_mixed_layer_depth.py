#!/usr/bin/env python
"""
Mixing down to Mixed Layer Depth
================================
"""

from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers.reader_constant import Reader as ConstantReader

#%%
# Mixed Layer Depth of 20m West of 3 deg E, and 50m to the east
r1 = ConstantReader({'ocean_mixed_layer_thickness': 20})
r2 = ConstantReader({'ocean_mixed_layer_thickness': 50})
r1.xmax = 3
r2.xmin = 3

#%%
# First with Sundby1983 parameterization of diffusivity, based on wind and MLD
o = OceanDrift(loglevel=50)
o.add_reader([r1, r2])
o.set_config('environment:constant:y_wind', 8)  # Some wind for mixing
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983')
# Increasing background diffusivity beyond default (1.2e-5) to avoid artefact due to sharp gradient at MLD
o.set_config('vertical_mixing:background_diffusivity', 0.001)
o.seed_cone(lon=[2, 4], lat=[60, 60], time=datetime.now(), number=5000)
o.run(duration=timedelta(hours=48))
o.animation_profile()

#%%
# .. image:: /gallery/animations/example_mixed_layer_depth_0.gif

#%%
# Same, but with Large1994 parameterization of diffusivity
o = OceanDrift(loglevel=50)
o.add_reader([r1, r2])
o.set_config('environment:constant:y_wind', 8)  # Some wind for mixing
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Large1994')
o.set_config('vertical_mixing:background_diffusivity', 0.001)
o.seed_cone(lon=[2, 4], lat=[60, 60], time=datetime.now(), number=5000)
o.run(duration=timedelta(hours=48))
o.animation_profile()

#%%
# .. image:: /gallery/animations/example_mixed_layer_depth_1.gif

#%%
# Using Large1994, but with 0 diffusivity below MLD
o = OceanDrift(loglevel=50)
o.add_reader([r1, r2])
o.set_config('environment:constant:y_wind', 8)  # Some wind for mixing
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Large1994')
o.set_config('vertical_mixing:background_diffusivity', 0)
o.seed_cone(lon=[2, 4], lat=[60, 60], time=datetime.now(), number=5000)
o.run(duration=timedelta(hours=48))
o.animation_profile()

#%%
# .. image:: /gallery/animations/example_mixed_layer_depth_2.gif
