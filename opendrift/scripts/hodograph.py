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
import matplotlib.pyplot as plt

sys.path.append("..")

try:
    from opendrift.readers import reader_netCDF_CF_generic
    from opendrift.readers import reader_ROMS_native
    readers = [reader_netCDF_CF_generic, reader_ROMS_native]
except ImportError: # development
    sys.exit('Please add opendrift folder to your PYTHONPATH.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',
                        help='<URL or netCDF filename>')
    parser.add_argument('lon', help='longitude in degrees')
    parser.add_argument('lat', help='longitude in degrees')
    parser.add_argument('time', help='YYYYMMDDHHMM')
    parser.add_argument('-e', action='store_true',
                        help='Report errors on failure.')

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

    time = datetime.strptime(args.time, '%Y%m%d%H%M')
    lon = np.atleast_1d(args.lon)
    lat = np.atleast_1d(args.lat)

    x, y = r.lonlat2xy(lon, lat)
    r.buffer = 3
    i=3; j=3  # center of block

    uv = r.get_variables(['x_sea_water_velocity', 'y_sea_water_velocity',
                         'sea_floor_depth_below_sea_level'],
                         time, x, y, r.z)

    u_comp = uv['x_sea_water_velocity'][:,i,j]
    v_comp = uv['y_sea_water_velocity'][:,i,j]
    depth = uv['sea_floor_depth_below_sea_level'][i,j]
    u_comp = u_comp[0:sum(~u_comp.mask)]
    v_comp = v_comp[0:sum(~v_comp.mask)]

    u_rot, v_rot = r.rotate_vectors(x, y, u_comp, v_comp,
                                    r.proj, '+proj=latlong')

    vel = np.sqrt(u_rot**2 + v_rot**2)
    mx = np.nanmax((np.abs(u_rot), np.abs(v_rot)))
    fig = plt.figure()
    ax = fig.gca()
    plt.quiver(np.zeros(len(u_rot)), np.zeros(len(u_rot)), u_rot, v_rot, r.z[0:sum(~u_rot.mask)],
               scale=1.1, scale_units='x')
    plt.axis([-mx, mx, -mx, mx])
    plt.grid('on')
    cb = plt.colorbar(ticks=r.z[0:sum(~u_rot.mask)])
    cb.set_label('Depth [m]')
    plt.ylabel('North component [m/s]')
    plt.xlabel('East component [m/s]')
    textstr = u'%s UTC\n%.3fN\N{DEGREE SIGN}, %.3fE\N{DEGREE SIGN}' % \
             (time, float(args.lat), float(args.lon))
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    plt.show()

if __name__ == '__main__':
    main()

