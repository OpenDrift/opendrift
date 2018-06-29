#!/usr/bin/env python

# Example to illustrate stranding options using an artificial
# east-west oscillating current field
# Knut-Frode Dagestad, Feb 2017

from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_oscillating

from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Adding nordic reader for coastline
reader_nordic = reader_ROMS_native.Reader(o.test_data_folder() +
    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
reader_nordic.variables = ['land_binary_mask']

reader_osc = reader_oscillating.Reader('x_sea_water_velocity', amplitude=1,
                                       zero_time=reader_nordic.start_time)
o.add_reader([reader_osc])  # Oscillating east-west current component

o.fallback_values['y_sea_water_velocity'] = .2  # Adding northwards drift
o.set_config('general:basemap_resolution', 'i')

##########################################################
# Try different options: 'previous', 'stranding', 'none'
o.set_config('general:coastline_action', 'previous')
##########################################################

time = reader_osc.zero_time
lon = 12.2; lat = 67.7
o.seed_elements(lon, lat, radius=5000, number=15, time=time)

o.run(steps=36*4, time_step=900)

# Print and plot results
print(o)
o.plot()
#o.plot(background=['x_sea_water_velocity', 'y_sea_water_velocity'])
o.animation()
