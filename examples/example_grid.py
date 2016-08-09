#!/usr/bin/env python

from datetime import datetime, timedelta

import numpy as np

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil
from models.oceandrift import OceanDrift

#o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information
o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('../test_data/16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst = reader_netCDF_CF_generic.Reader('../test_data/16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

# GlobCurrent
#reader_globcurrent = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREUL_HS-ALT_SUM-V01.0_FULL_TIME_SERIE')  # Total
#reader_globcurrent = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/GC_MOD_TIDE_GLO_010_FES2012_FULL_TIME_SERIE') # FES Tidal

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=3.5, llcrnrlat=59.9,
                    urcrnrlon=5.5, urcrnrlat=61.2,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_norkyst])

# Seeding some particles
lons = np.linspace(3.5, 5.0, 100)
lats = np.linspace(60, 61, 100)
lons, lats = np.meshgrid(lons, lats)
lons = lons.ravel()
lats = lats.ravel()

# Seed oil elements at defined position and time
o.seed_elements(lons, lats, radius=0, number=10000,
                time=reader_norkyst.start_time)

print o

# Running model (until end of driver data)
o.run(steps=66*4, time_step=900)

# Print and plot results
print o
o.animation()
o.plot()
