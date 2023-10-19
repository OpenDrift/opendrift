#!/usr/bin/env python
"""
Seafloor interaction
=====================
"""

#%%
# This example demonstrates the three possibilities for
# interaction of particles with seafloor:
# - 'previous': particles are moved back to previous location
# - 'deactivate': particles are deactivated
# - 'lift_to_seafloor': particles are lifted vertically to seafloor level
#
# This is controlled by the config setting:
# o.set_config('general:seafloor_action', <action>)


from datetime import timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_oscillating


# readers
o = OceanDrift(loglevel=50)
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
reader_osc = reader_oscillating.Reader('x_sea_water_velocity', amplitude=10, period_seconds=3600)

runs = []
seafloor_actions = ['previous', 'deactivate', 'lift_to_seafloor']

for seafloor_action in seafloor_actions:

    o = OceanDrift(loglevel=50)  # Set loglevel to 0 for debug information
    o.set_config('drift:max_speed', 10)

    o.add_reader([reader_osc, reader_norkyst])

    o.set_config('drift:horizontal_diffusivity', 0)
    o.set_config('environment:constant:y_sea_water_velocity', 0)
    o.set_config('environment:constant:land_binary_mask', 0)
    o.set_config('general:use_auto_landmask', False)
    o.set_config('general:seafloor_action', seafloor_action)

    # Seeding some particles 50m above seafloor
    o.seed_elements(lon=4.2, lat=62.0, z='seafloor+50',
                    radius=7000, number=200, time=reader_norkyst.start_time)

    o.run(duration=timedelta(hours=2), time_step=120)
    runs.append(o)


runs[0].animation(fast=False, compare=runs[1:], legend=seafloor_actions,
                  vmin=0, vmax=300, background='sea_floor_depth_below_sea_level')
#%%
# .. image:: /gallery/animations/example_seafloor_interaction_0.gif

#%%
# For the run with 'lift_to_seafloor', we see that
# elements have been lifted vertically
runs[2].plot_property('z')
