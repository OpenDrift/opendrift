#!/usr/bin/env python

from datetime import timedelta

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil3D import OpenOil3D

o = OpenOil3D(loglevel=0)  # Set loglevel to 0 for debug information

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

o.add_reader([reader_norkyst])
o.fallback_values['x_wind'] = 0
o.fallback_values['y_wind'] = 0

# Seeding some particles
lon = 4.5; lat = 62.0
time = [reader_norkyst.start_time,
        reader_norkyst.start_time + timedelta(hours=1)]

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, z='seafloor', radius=5, number=3000, time=time)

o.set_config('processes:turbulentmixing', True)  # Otherwise also no updrift

# Running model
o.run(steps=12*2, time_step=300, time_step_output=300)

# Print and plot results
print(o)
o.animation_profile()
