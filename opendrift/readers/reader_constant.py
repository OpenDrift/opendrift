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

from opendrift.readers.basereader import BaseReader, ContinuousReader
import numpy as np


class Reader(BaseReader, ContinuousReader):
    '''A very simple reader that always give the same value for its variables'''

    def __init__(self, parameter_value_map):
        '''init with a map {'variable_name': value, ...}'''

        for key, var in parameter_value_map.items():
            parameter_value_map[key] = np.atleast_1d(var)
        self._parameter_value_map = parameter_value_map
        self.variables = list(parameter_value_map.keys())
        self.proj4 = '+proj=latlong'
        self.xmin = -180
        self.xmax = 180
        self.ymin = -90
        self.ymax = 90
        self.start_time = None
        self.end_time = None
        self.time_step = None
        self.name = 'constant_reader'

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, z=None):

        variables = {'time': time, 'x': x, 'y': y, 'z': z}
        #variables.update(self._parameter_value_map)
        for var in requestedVariables:
            variables[var] = self._parameter_value_map[var]*np.ones(x.shape)

        return variables

