#!/usr/bin/env python

# Seeding elements around the border of a ocean model domain (NorKyst800)
# to demonstrate autmatic transition back and forth with another model
# covering a larger domain (Nordic)

import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

# Nordic4
reader_nordic4 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/nordic4km/zdepths1h/aggregate_be')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=9.5, llcrnrlat=68.8,
                    urcrnrlon=19.0, urcrnrlat=71.2,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_norkyst, reader_nordic4])
#o.add_reader([reader_basemap, reader_norkyst])

# Seeding some particles
lons = np.linspace(10.2, 12.2, 50)
lats = np.linspace(69.8, 70.8, 50)
lons, lats = np.meshgrid(lons, lats)
lons = lons.ravel()
lats = lats.ravel()

# Seed oil elements at defined position and time
o.seed_elements(lons, lats, radius=0, number=2500,
                time=reader_nordic4.start_time)

# Running model (until end of driver data)
o.run(steps=16*4, time_step=900)

# Print and plot results
print(o)
o.animation()
