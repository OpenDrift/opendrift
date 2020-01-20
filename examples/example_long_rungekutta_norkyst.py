#!/usr/bin/env python
"""
Runge-Kutta scheme on Norkyst model
==================================

Illustrating the difference between Euler and Runge-Kutta propagation
schemes, using a "real" current fields from the NorKyst800 model
"""

from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
time = reader_norkyst.start_time
reader_norkyst.interpolation = 'linearND'

reader_landmask = reader_global_landmask.Reader(
                    llcrnrlon=4, llcrnrlat=59.9,
                    urcrnrlon=5.5, urcrnrlat=61.5, resolution='h')

o.add_reader([reader_norkyst, reader_landmask])
lon = 4.5; lat = 60.0;

# First run, with Euler scheme:
o.set_config('drift:scheme', 'euler')
o.seed_elements(lon, lat, radius=0, number=1, time=time)
o.run(steps=66*2, time_step=1800)

# Second run, with Runge-Kutta scheme:
o2 = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
o2.add_reader([reader_norkyst, reader_landmask])
o2.set_config('drift:scheme', 'runge-kutta')
o2.seed_elements(lon, lat, radius=0, number=1, time=time)
o2.run(steps=66*2, time_step=1800)

# Animate and compare the two runs
o.animation(compare=o2, legend=['Euler scheme', 'Runge-Kutta scheme'], filename='rungekutta_norkyst.gif')

#%%
# .. image:: /gallery/animations/rungekutta_norkyst.gif
