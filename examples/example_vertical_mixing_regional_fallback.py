#!/usr/bin/env python
"""
Vertical mixing with fallback diffusivity model
===============================================

With ``diffusivitymodel='environment'`` and ``diffusivitymodel_fallback``,
elements use the reader diffusivity where it has coverage, and the fallback
model elsewhere.
"""

from datetime import datetime, timedelta
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers.reader_constant import Reader as ConstantReader

#%%
# Diffusivity from an ocean model west of 3°E only, and an eastward current
r_ocean_model = ConstantReader({'ocean_vertical_diffusivity': 0.01})
r_ocean_model.xmax = 3
r_current = ConstantReader({'x_sea_water_velocity': 0.3, 'y_sea_water_velocity': 0})

o = OceanDrift(loglevel=20)
o.add_reader([r_ocean_model, r_current])
o.set_config('environment:constant:y_wind', 8)
o.set_config('environment:fallback:ocean_mixed_layer_thickness', 30)
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'environment')
o.set_config('vertical_mixing:diffusivitymodel_fallback', 'windspeed_Large1994')
o.set_config('vertical_mixing:background_diffusivity', 0.0001)
o.seed_cone(lon=[2, 2.5], lat=[60, 60], time=datetime.now(), number=2000, z=-5)
o.run(duration=timedelta(hours=72), time_step=1800)
o.animation_profile()

#%%
# .. image:: /gallery/animations/example_vertical_mixing_regional_fallback_0.gif
#
# Elements east of 3°E are mixed with the Large1994 diffusivity (wind and
# mixed layer depth), the others with the ocean model diffusivity.
