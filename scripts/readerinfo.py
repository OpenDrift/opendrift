#!/usr/bin/env python

# Utility script to display contents of a netCDF CF-compliant
# file or URL containing driver data suitable for opendrift
#
# Knut-Frode Dagestad, 19 Feb 2015

import sys
import argparse

try:
    from readers import reader_netCDF_CF_generic
except ImportError: # development
    sys.exit('Please add opendrift folder to your PYTHONPATH.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',
                        help='<URL or netCDF filename>')
    parser.add_argument('-p', dest='variable',
                        default='noplot', nargs='?',
                        help='Plot domain (or variable if given)')

    args = parser.parse_args()

    r = reader_netCDF_CF_generic.Reader(args.filename)
    print r

    if args.variable != 'noplot':
        if args.variable is None:
            r.plot()
        else:
            r.plot(args.variable)
