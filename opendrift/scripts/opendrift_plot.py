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

import sys
import argparse
import numpy as np
sys.path.append("..")

import opendrift

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',
                        help='<OpenDrift output filename (netCDF)>')
    parser.add_argument('-b', dest='buffer',
                        default=0.1,
                        help='Buffer around plot in degrees lon/lat.')
    parser.add_argument('-f', dest='outfile',
                        default=None, help='Save plot to filename.')

    args = parser.parse_args()


    o = opendrift.open(args.filename)
    print(o)
    o.plot(buffer=float(args.buffer), filename=args.outfile)

if __name__ == '__main__':
    main()
