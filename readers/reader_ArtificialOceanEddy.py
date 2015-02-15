#!/usr/bin/env python

import os
from datetime import datetime, timedelta

import numpy as np
from netCDF4 import Dataset

from readers import Reader, pyproj


class Reader(Reader):

    def __init__(self, lon=2, lat=66, proj4='+proj=stere +lat_0=90 +lon_0=70 +lat_ts=60 +units=m +a=6.371e+06 +e=0 +no_defs'):
    
        self.fileName = 'ArtificialOceanEddy'
        self.name = 'ArtificialOceanEddy'

        self.proj4 = proj4
        self.proj = pyproj.Proj(proj4)

        # Calculate x,y of center of eddy from given lon and lat
        x0, y0 = self.proj(lon, lat)
        width = 300000  # 200 km box
        print x0, y0
        self.xmin = x0 - width
        self.xmax = x0 + width
        self.ymin = y0 - width
        self.ymax = y0 + width
        self.delta_x = 1000  # Artificial 1 km pixel size
        self.delta_y = 1000  # Artificial 1 km pixel size
        self.startTime = None
        self.endTime = None
        self.timeStep = None
        self.variables = ['x_sea_water_velocity', 'y_sea_water_velocity']


        # Run constructor of parent Reader class
        super(Reader, self).__init__()


    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, depth=None, block=False):

        requestedVariables, time, x, y, depth, outside = self.check_arguments(
            requestedVariables, time, x, y, depth)

        # TBD

        return variables
