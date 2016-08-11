#!/usr/bin/env python

import os
from datetime import datetime, timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

if not os.path.exists('openoil.nc'):
    raise ValueError('Please run example.py first to generate a '
                     'netCDF file to be imported.')

o.io_import_file('openoil.nc')

print o

o.plot(buffer=.1)
o.plot_property('mass_oil')
