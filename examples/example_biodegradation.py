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
o.set_config('drift', {'current_uncertainty': 0, 'wind_uncertainty': 0})

#%%
# Seeding some particles
o.set_config('drift:vertical_mixing', False)  # Keeping oil at fixed depths
o.set_config('processes:biodegradation', True)
o.set_config('biodegradation:method', 'decay_rate')

#%%
# Fast decay for droplets, and slow decay for slick 
decay = {'biodegradation_decay_rate_slick': np.log(2)/timedelta(days=100).total_seconds(),
         'biodegradation_decay_rate_droplet': np.log(2)/timedelta(days=3).total_seconds(),
         'oil_type': 'GENERIC MEDIUM CRUDE', 'm3_per_hour': .5}

#%%
# Seed 5000 oil elements at surface, and 5000 elements at 100m depth
o.seed_elements(lon=4, lat=60.0, z=0, number=5000, time=datetime.now(), **decay)
o.seed_elements(lon=4, lat=60.0, z=-100, number=5000, time=datetime.now(), **decay)

#%%
# Running model
o.run(duration=timedelta(hours=24), time_step=900) 

#%%
# Print and plot results
o.plot_oil_budget(show_watercontent_and_viscosity=False, show_wind_and_current=False)
