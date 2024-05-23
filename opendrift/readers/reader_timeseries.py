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
# Copyright 2020, Knut-Frode Dagestad, MET Norway

from opendrift.readers.basereader import BaseReader, ContinuousReader
import numpy as np

class Reader(BaseReader, ContinuousReader):
    '''Reader providing the nearest value in time from a given time series.'''

    def __init__(self, parameter_value_map):
        '''init with a map {'time':, time array, 'variable_name': value, ...}
        If there is the key lon or lat in the map, it will be stored as
        self.lon and self.lat but not as a timeserie.
        '''

        self.times = parameter_value_map['time']
        if type(self.times) is not list:
            raise ValueError('time must be a list of datetime objects')
        del parameter_value_map['time']

        for key, var in parameter_value_map.items():
            if key == 'lon':
                self.lon=var
                continue
            if key == 'lat':
                self.lat=var
                continue
            parameter_value_map[key] = np.atleast_1d(var)
            if len(parameter_value_map[key]) != len(self.times):
                raise ValueError('All variables must have same length as time array')

        self._parameter_value_map = parameter_value_map
        self.variables = list(parameter_value_map.keys())
        self.proj4 = '+proj=latlong'
        self.xmin = -180
        self.xmax = 180
        self.ymin = -90
        self.ymax = 90
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        #self.end_time = None
        self.time_step = self.times[1] - self.times[0]
        self.name = 'reader_timeseries'

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, z=None):

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        variables = {'time': time, 'x': x, 'y': y, 'z': z}
        #variables.update(self._parameter_value_map)
        for var in requestedVariables:
            variables[var] = self._parameter_value_map[var][indxTime]*np.ones(x.shape)

        return variables

