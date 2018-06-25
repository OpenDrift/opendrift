#!/usr/bin/env python

import os

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic

from opendrift.models.oceandrift import OceanDrift

ncfile = 'backandforth.nc'

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=3, llcrnrlat=59,
                urcrnrlon=6, urcrnrlat=62, resolution='h')

o.add_reader([reader_norkyst, reader_basemap])

################
# Forward run
################
# Seeding some particles
lon = 4.2; lat = 60.1; 
time = reader_norkyst.start_time
o.seed_elements(lon, lat, radius=1000, number=100, time=time)

o.run(steps=50*4, time_step=900, outfile=ncfile)

# Print and plot results
print(o)
o.plot(buffer=.2, filename='forward.png')

#################
## Backward run
#################
o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
#o.add_reader([reader_norkyst, reader_basemap])
# Import forward run, and seed elements at final positions and time
o.io_import_file(ncfile)
elements_final = o.elements
time_final = o.time
del o
o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
o.add_reader([reader_norkyst, reader_basemap])
o.schedule_elements(elements_final, time_final)

# Running model backwards from end of forward simulation
o.run(steps=50*4, time_step=-900)

# Print and plot results
print(o)
o.plot(buffer=.2, filename='backward.png')
os.remove(ncfile)
print('#############################################')
print('Compare plots forward.png and backward.png')
print('#############################################')
