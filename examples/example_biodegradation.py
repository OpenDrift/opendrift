#!/usr/bin/env python
"""
Biodegradation of oil
======================
"""

import numpy as np
from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information
time = datetime.now()

# No motion is needed for this test
o.set_config('environment:constant', {k: 0 for k in
             ['x_wind', 'y_wind', 'x_sea_water_velocity', 'y_sea_water_velocity']})
o.set_config('drift', {'current_uncertainty': 0, 'wind_uncertainty': 0, 'horizontal_diffusivity': 10})

#%%
# Seeding some particles
o.set_config('drift:vertical_mixing', True)
o.set_config('processes:biodegradation', True)
o.set_config('biodegradation:method', 'half_time')

#%%
# Fast decay for droplets, and slow decay for slick 
decay = {'biodegradation_half_time_slick': 5, # days
         'biodegradation_half_time_droplet': 1, # days,
         'oil_type': 'GENERIC MEDIUM CRUDE', 'm3_per_hour': .5, 'diameter': 1e-5}  # small droplets

#%%
# Seed 500 oil elements at surface, and 500 elements at 50m depth
o.seed_elements(lon=4, lat=60.0, z=0, number=500, time=datetime.now(), **decay)
o.seed_elements(lon=4, lat=60.0, z=-50, number=500, time=datetime.now(), **decay)

#%%
# Running model
o.run(duration=timedelta(hours=72), time_step=3600)

#%%
# Plot results
o.plot_oil_budget(show_watercontent_and_viscosity=False, show_wind_and_current=False)

o.animation_profile(markersize='mass_oil', markersize_scaling=50, color='z', alpha=.5)

#%%
# .. image:: /gallery/animations/example_biodegradation_0.gif
