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

from datetime import datetime, timedelta
import numpy as np
import pyproj

from opendrift.readers.basereader import BaseReader, ContinuousReader


class Reader(BaseReader, ContinuousReader):
    """Analytical current field with double gyre"""

    def __init__(self, initial_time=datetime(2000,1,1,0,0),
                 epsilon=0.1, omega=0.628, A=0.25,
                 proj4='+proj=stere +lat_0=0 +lon_0=0 +lat_ts=0 ' +
                       '+units=m +a=6.371e+06 +e=0 +no_defs'):

        self.fileName = 'double_gyre'
        self.name = 'double_gyre'

        self.proj4 = proj4
        self.proj = pyproj.Proj(proj4)

        self.xmin = 0.
        self.xmax = 2.
        self.ymin = 0.
        self.ymax = 1.
        self.A = A
        self.epsilon=epsilon
        self.omega=omega
        self.initial_time = initial_time
        self.start_time = None
        self.end_time = None
        self.time_step = None
        self.variables = ['x_sea_water_velocity',
                          'y_sea_water_velocity']

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, z=None):

        requestedVariables, time, x, y, z, outside = self.check_arguments(
            requestedVariables, time, x, y, z)

        variables = {}

        # Construct double gyre current field
        # https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html
        t = (time-self.initial_time).total_seconds()
        a = self.epsilon*np.sin(self.omega*t)
        b = 1-2*self.epsilon*np.sin(self.omega*t)
        f = a*x*x + b*x
        dfdx = 2*a*x + b

        variables['x_sea_water_velocity'] = -np.pi*self.A*np.sin(np.pi*f)*np.cos(np.pi*y)
        variables['y_sea_water_velocity'] = np.pi*self.A*np.cos(np.pi*f)*np.sin(np.pi*y)*dfdx

        variables['land_binary_mask'] = np.zeros(x.shape)
        variables['time'] = time
        variables['x'] = x
        variables['y'] = y

        return variables
