#!/usr/bin/env python

# Utility script to display contents of a netCDF CF-compliant
# file or URL containing driver data suitable for opendrift
#
# Knut-Frode Dagestad, 19 Feb 2015

import sys

try:
    from readers import reader_netCDF_CF_generic
except ImportError: # development
    sys.exit('Please add opendrift folder to your PYTHONPATH.')

if (len(sys.argv) != 2):
    sys.exit('Usage: readerinfo <filename or URL>')

r = reader_netCDF_CF_generic.Reader(sys.argv[1])
print r
