#!/usr/bin/env python

import argparse
import numpy as np
from models.openoil import OpenOil

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',
                        help='<OpenDriftaoutput filename (netCDF)>')
    parser.add_argument('-b', dest='buffer',
                        default=1.0,
                        help='Buffer around plot in degrees lon/lat.')

    args = parser.parse_args()


    o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

    print o
    print type(np.array(args.buffer))
    o.io_import_file(args.filename)
    o.plot(buffer=np.float(args.buffer))
