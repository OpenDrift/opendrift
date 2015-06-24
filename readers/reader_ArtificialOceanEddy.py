#!/usr/bin/env python

from datetime import datetime, timedelta

import numpy as np

from readers import Reader, pyproj


class Reader(Reader):
    """Artificial reader, with cyclonic surface current around selected centre.

    Purpose is demonstration and testing (unittest).
    """

    def __init__(self, lon=2, lat=66,
                 proj4='+proj=stere +lat_0=90 +lon_0=70 +lat_ts=60 ' +
                       '+units=m +a=6.371e+06 +e=0 +no_defs'):

        self.fileName = 'ArtificialOceanEddy'
        self.name = 'ArtificialOceanEddy'

        self.proj4 = proj4
        self.proj = pyproj.Proj(proj4)

        # Calculate x,y of center of eddy from given lon and lat
        x0, y0 = self.proj(lon, lat)
        width = 600000  # 600 km box
        self.pixelsize = 10000  # Artificial 10 km pixel size
        self.delta_x = self.pixelsize
        self.delta_y = self.pixelsize
        self.x0 = x0
        self.y0 = y0
        self.xmin = x0 - width
        self.xmax = x0 + width
        self.ymin = y0 - width
        self.ymax = y0 + width
        self.start_time = None
        self.end_time = None
        self.time_step = None
        self.variables = ['x_sea_water_velocity', 'y_sea_water_velocity']

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, depth=None, block=False):

        requestedVariables, time, x, y, depth, outside = self.check_arguments(
            requestedVariables, time, x, y, depth)

        variables = {}
        # Construct cyclonic current field
        size = np.int(self.xmax-self.xmin) / self.pixelsize + 1
        if block is True:
            x = np.linspace(self.xmin, self.xmax, size)
            y = np.linspace(self.ymin, self.ymax, size)
            X, Y = np.meshgrid(x-self.x0, y-self.y0)
        else:
            X = x - self.x0
            Y = y - self.y0
        radius = np.sqrt(X*X + Y*Y)
        radius[radius == 0] = 1

        variables['time'] = time
        variables['x'] = x
        variables['y'] = y
        variables['x_sea_water_velocity'] = -Y / radius  # np.ones((size,size))
        variables['y_sea_water_velocity'] = X / radius  # np.zeros((size,size))

        return variables
