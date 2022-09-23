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

import argparse
import numpy as np
import opendrift

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename1',
                        help='<OpenDrift filename (netCDF) of first file>')
    parser.add_argument(dest='filename2', nargs='?',
                        help='<OpenDrift filename (netCDF) of second file>')
    parser.add_argument('-b', dest='buffer',
                        default=1.0,
                        help='Buffer around plot in degrees lon/lat.')
    parser.add_argument('-f', dest='outfile',
                        default=None, help='Save animation to filename.')
    parser.add_argument('-c', dest='color',
                        default=False, help='Color elements with this parameter.')

    args = parser.parse_args()


    o1 = opendrift.open(args.filename1)
    print(o1)

    if args.filename2 is None:
        o1.animation(filename=args.outfile, color=args.color)
    else:
        o2 = opendrift.open(args.filename2)
        print(o2)

        # Animate and compare the two runs
        o1.animation(compare=o2, legend=[args.filename1, args.filename2],
                     filename=args.outfile, color=args.color)

if __name__ == '__main__':
   main()
