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

import numpy as np
import pyproj

from opendrift.readers.basereader import BaseReader, ContinuousReader


class Reader(BaseReader, ContinuousReader):
    """Artificial reader, with cyclonic surface current around selected centre.

    Purpose is demonstration and testing (unittest).
    """

    def __init__(self, lon=2, lat=66,
                 proj4='+proj=stere +lat_0=90 +lon_0=70 +lat_ts=60 ' +
                       '+units=m +a=6.371e+06 +no_defs'):

        self.fileName = 'ArtificialOceanEddy'
        self.name = 'ArtificialOceanEddy'

        self.proj4 = proj4
        self.proj = pyproj.Proj(self.proj4)

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
                      x=None, y=None, z=None):

        requestedVariables, time, x, y, z, outside = self.check_arguments(
            requestedVariables, time, x, y, z)

        variables = {}
        # Construct cyclonic current field
        size = int(int(self.xmax-self.xmin) / self.pixelsize + 1)
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
