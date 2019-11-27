"""
Model landmask
===============================

Comparing two simulation runs, with and without wind
"""
#!/usr/bin/env python

# Comparing two simulation runs, with and without wind

from datetime import timedelta

from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift

lon = 14.75; lat = 68.1

o = OceanDrift(loglevel=0)

reader_nordic = reader_ROMS_native.Reader(o.test_data_folder() +
    '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')

# First run, with landmask
o.add_reader([reader_nordic])
time = reader_nordic.start_time
o.seed_elements(lon, lat, radius=3000, number=1000, time=time)
o.set_config('general:use_auto_landmask', True)
o.run(end_time=reader_nordic.end_time, time_step=1800)

# Second run, with landmask from ocean model
o2 = OceanDrift(loglevel=0)
o2.add_reader([reader_nordic])
lon = 14.75; lat = 68.1
o2.seed_elements(lon, lat, radius=3000, number=1000, time=time)
o2.set_config('general:use_auto_landmask', False)
o2.run(end_time=reader_nordic.end_time, time_step=1800)


# Animation illustrating that red particles strand at ocean model land cells, and black particles strand at Basemap land polygons
o.animation(compare=o2, background='land_binary_mask',
            legend=['Landmask', 'Ocean model landmask'])
