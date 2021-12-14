#!/usr/bin/env python
"""
Seafloor oil spill
===================
"""

from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader([reader_norkyst])
o.set_config('environment:fallback:x_wind', 3)
o.set_config('environment:fallback:y_wind', 7)
o.set_config('drift:vertical_mixing', True)

#%%
# Seeding some particles
time = [reader_norkyst.start_time,
        reader_norkyst.start_time + timedelta(hours=1)]
o.seed_elements(lon=4.5, lat=62.0, z='seafloor', radius=0, number=3000,
                time=time, oil_type='GENERIC DIESEL')

#%%
# Setting the range of droplet sizes for the seafloor release
o.set_config('seed:droplet_diameter_min_subsea', 0.0001)
o.set_config('seed:droplet_diameter_max_subsea', 0.0005)

#%%
# Running model with a small timestep to resolve the boyant rising
o.run(duration=timedelta(hours=2), time_step=60, time_step_output=60)

#%%
# Print and plot results
print(o)

o.animation_profile()
#%%
# .. image:: /gallery/animations/example_oilspill_seafloor_0.gif

o.animate_vertical_distribution(bins=30, subsamplingstep=5)
#%%
# .. image:: /gallery/animations/example_oilspill_seafloor_1.gif

o.plot_oil_budget()
