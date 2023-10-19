#!/usr/bin/env python
"""
Satellite
==================================
"""

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

#%%
# Adjusting some configuration
o.set_config('drift:vertical_mixing', False)
o.set_config('processes:dispersion', True)
o.set_config('processes:evaporation', False)
o.set_config('processes:emulsification', True)
o.set_config('drift:current_uncertainty', .1)  # Diffusion
o.set_config('drift:wind_uncertainty', 1)

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([reader_norkyst, reader_arome])

#%%
# Seed oil particles within contour detected from satellite
o.seed_from_gml(o.test_data_folder() + 'radarsat_oil_satellite_observation/RS2_20151116_002619_0127_SCNB_HH_SGF_433012_9730_12182143_Oil.gml',
    num_elements=2000)

#%%
# Running model for 6 hours
o.run(steps=6*4, time_step=900)

#%%
# Print and plot results
print(o)
o.animation(fast=True, buffer=0.1)

#%%
# .. image:: /gallery/animations/example_satellite_0.gif

o.plot(fast=True, buffer=0.1)
