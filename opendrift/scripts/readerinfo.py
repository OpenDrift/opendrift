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
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
# 
# Copyright 2015, Knut-Frode Dagestad, MET Norway


# Utility script to display contents of a netCDF CF-compliant
# file or URL containing driver data suitable for opendrift
#
# Knut-Frode Dagestad, 19 Feb 2015

import sys
import argparse
from datetime import datetime
import numpy as np

try:
    from opendrift.readers import reader_netCDF_CF_generic
    from opendrift.readers import reader_ROMS_native
    try:
        from opendrift.readers import reader_grib
        readers = [reader_netCDF_CF_generic, reader_ROMS_native,
                   reader_grib]
    except:
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
    parser.add_argument('-vmin', dest='vmin',
                        default=None, nargs='?',
                        help='Minimum value for colorbar')
    parser.add_argument('-vmax', dest='vmax',
                        default=None, nargs='?',
                        help='Maximum value for colorbar')
    parser.add_argument('-e', action='store_true',
                        help='Report errors on failure.')
    parser.add_argument('-lon', dest='lon', default=None,
                        help='Report data from position.')
    parser.add_argument('-lat', dest='lat', default=None,
                        help='Report data from position.')
    parser.add_argument('-time', dest='time', default=None,
                        help='Report data from position at time [YYYYmmddHHMM].')

    args = parser.parse_args()

    for reader in readers:
        try:
            print('Testing %s...' % reader.__file__)
            r = reader.Reader(args.filename)
            print(r)
            break
        except Exception as me:
            if args.e is True:
                print(me)
                import traceback
                print(traceback.format_exc())
                print('---------------------------------------')
            print('...not applicable.')

    if not 'r' in locals():            
        sys.exit('No readers applicable for ' + args.filename)

    if args.lon is not None:
        if args.lon is None or args.time is None:
            raise ValueError('Both lon, lat and time must be given')
        lon = np.atleast_1d(float(args.lon))
        lat = np.atleast_1d(float(args.lat))
        x,y = r.lonlat2xy(lon, lat)
        time = datetime.strptime(args.time, '%Y%m%d%H%M')
        r.buffer=3
        i=3; j=3  # center of block
        variables = [var for var in r.variables if
                     var not in ('time') and 'time' not in var]
        data = r.get_variables(variables, time, x, y, z=0, block=True)
        for var in variables:
            print('%s : %s' % (var, data[var][i,j]))
        if 'x_wind' in variables and 'y_wind' in variables:
            print('windspeed : %s' % np.sqrt(
                data['x_wind'][i,j]**2 + data['y_wind'][i,j]**2))

    if args.variable != 'noplot':
        if args.variable is None:
            r.plot()
        else:
            if args.vmin is None:
                vmin = None
            else:
                vmin = np.float(args.vmin)
            if args.vmax is None:
                vmax = None
            else:
                vmax = np.float(args.vmax)
            
            r.plot(args.variable, vmin=vmin, vmax=vmax)
