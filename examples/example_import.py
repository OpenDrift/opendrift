#!/usr/bin/env python

import os

import opendrift


if not os.path.exists('openoil.nc'):
    raise ValueError('Please run example.py first to generate a '
                     'netCDF file to be imported.')

o = opendrift.open('openoil.nc')
print(o)

o.plot()
o.plot_property('mass_oil')
