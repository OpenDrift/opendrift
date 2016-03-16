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

import argparse
import numpy as np
from models.openoil3D import OpenOil3D

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename1',
                        help='<OpenDrift filename (netCDF) of first file>')
    parser.add_argument(dest='filename2', nargs='?',
                        help='<OpenDrift filename (netCDF) of second file>')
    parser.add_argument('-b', dest='buffer',
                        default=1.0,
                        help='Buffer around plot in degrees lon/lat.')

    args = parser.parse_args()


    o1 = OpenOil3D(loglevel=0)  # Set loglevel to 0 for debug information
    o1.io_import_file(args.filename1)
    print o1

    if args.filename2 is None:
        o1.animation()
    else:
        o2 = OpenOil3D(loglevel=0)  # Set loglevel to 0 for debug information
        o2.io_import_file(args.filename2)
        print o2

        # Animate and compare the two runs
        o1.animation(compare=o2, legend=[args.filename1, args.filename2])
