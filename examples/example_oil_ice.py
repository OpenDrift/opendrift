#!/usr/bin/env python
"""
Oil in ice
==================================
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil
import cmocean

o = OpenOil(loglevel=20)

#%%
# Using live data from Barents 2.5 km ocean model
o.add_readers_from_list(['https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'])

#%%
# Adjusting some configuration
o.set_config('processes:dispersion',  False)
o.set_config('processes:evaporation',  False)
o.set_config('processes:emulsification',  False)
o.set_config('drift:horizontal_diffusivity', 10)
o.set_config('drift:truncate_ocean_model_below_m', 3)

#%%
# Imaginary oil spill in Hinlopen strait
o.seed_elements(lon=19.1909, lat=79.5986, radius=50,
                number=3000, time=datetime.utcnow() - timedelta(days=7))

#%%
# Running model
o.run(duration=timedelta(hours=48), time_step=1800, time_step_output=3600)

#%%
# Print and plot results
print(o)
o.animation(background='sea_ice_area_fraction', cmap=cmocean.cm.ice,
            vmin=0, vmax=1, bgalpha=1, fast=False)

#%%
# .. image:: /gallery/animations/example_oil_ice_0.gif

o.plot(background='sea_ice_area_fraction', cmap=cmocean.cm.ice,
       vmin=0, vmax=1, bgalpha=1, fast=False)
