#!/usr/bin/env python
#
# This file is part of OpenDrift.
# 
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
# 
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
# 
# Copyright 2015, Knut-Frode Dagestad, MET Norway


# Utility script to display contents of a netCDF CF-compliant
# file or URL containing driver data suitable for opendrift
#
# Knut-Frode Dagestad, 19 Feb 2015

import sys
import argparse

try:
    from readers import reader_netCDF_CF_generic
    from readers import reader_ROMS_native
    readers = [reader_netCDF_CF_generic, reader_ROMS_native]
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

    for reader in readers:
        try:
            r = reader.Reader(args.filename)
            print r
            break
        except Exception as me:
            print me
            import traceback
            print traceback.format_exc()
            print 'Reader not applicable: ' + str(reader)
            print '---------------------------------------'

    if not 'r' in locals():            
        sys.exit('No readers applicable for ' + args.filename)

    if args.variable != 'noplot':
        if args.variable is None:
            r.plot()
        else:
            r.plot(args.variable)
