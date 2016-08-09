#!/usr/bin/env python

# Illustrating the difference between Euler and Runge-Kutta propagation
# schemes, using a "real" current fields from the NorKyst800 model

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

reader_norkyst = reader_netCDF_CF_generic.Reader(
    '../test_data/16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
time = reader_norkyst.start_time
reader_norkyst.interpolation = 'linearND'

reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=4, llcrnrlat=59.9,
                    urcrnrlon=5.5, urcrnrlat=61.5, resolution='h')

o.add_reader([reader_norkyst, reader_basemap])
lon = 4.5; lat = 60.0;

# First run, with Euler scheme:
o.config['drift']['scheme'] = 'euler'
o.seed_elements(lon, lat, radius=0, number=1, time=time)
o.run(steps=66*2, time_step=1800)

# Second run, with Runge-Kutta scheme:
o2 = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o2.add_reader([reader_norkyst, reader_basemap])
o2.config['drift']['scheme'] = 'runge-kutta'
o2.seed_elements(lon, lat, radius=0, number=1, time=time)
o2.run(steps=66*2, time_step=1800)

# Animate and compare the two runs
o.animation(compare=o2, legend=['Euler scheme', 'Runge-Kutta scheme'],
            markersize=10)
