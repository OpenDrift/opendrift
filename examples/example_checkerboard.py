#!/usr/bin/env python

import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('test_data/16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=3.5, llcrnrlat=59.9,
                    urcrnrlon=5.5, urcrnrlat=61.2,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_norkyst])

# Seeding particles in a checkerboard pattern
di = 5 # Horizontal number of particles per square 
dj = 5 # Vertical number of particles per square
lons = np.linspace(3.5, 5.0, 100)
lats = np.linspace(60, 61, 100)

ii = np.arange(len(lons))//di
jj = np.arange(len(lats))//dj
ii, jj = np.meshgrid(ii, jj)
board = (ii+jj)%2 > 0

lons, lats = np.meshgrid(lons, lats)

lons = lons[board].ravel()
lats = lats[board].ravel()

# Seed oil elements at defined position and time
o.seed_elements(lons, lats, radius=0, number=10000,
                time=reader_norkyst.start_time)

# Running model (until end of driver data)
o.run(steps=66*2, time_step=1800)

# Print and plot results
print(o)
o.animation(filename="example_checkerboard_anim.gif")
#o.plot()
