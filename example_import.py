#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

o.io_import_file('openoil.nc')

print o

o.plot(buffer=.1)
o.plot_property('mass_oil')
