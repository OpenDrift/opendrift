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

    def __init__(self, *args, **kwargs):
        """init with a map {'variable_name': value, ...}
        
        or as variable_name1=value1, variable_name2=value2, ...
        value can also be an array, and in this case the map/dictionary
        must also include `element_ID` which corresponds to the elements that
        shall receive the actual value:
            self.environment.<variable_name> --> value[element_ID = self.elements.ID]  (pseudo code)
        This is however more simply achived by specifying environment when seeding, see:
        https://opendrift.github.io/gallery/example_element_dependent_environment.html

        """

        if len(args) == 1 and isinstance(args[0], dict):
            parameter_value_map = args[0]
        elif kwargs:
            parameter_value_map = kwargs
        for key, var in parameter_value_map.items():
            parameter_value_map[key] = np.atleast_1d(var)
        self._parameter_value_map = parameter_value_map
        self.variables = list([v for v in parameter_value_map.keys() if v !='element_ID'])
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

        if 'element_ID' in parameter_value_map:
            self._element_ID = True  # will be updated with indices of actual elements

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, z=None):

        variables = {'time': time, 'x': x, 'y': y, 'z': z}

        for var in requestedVariables:
            variables[var] = np.nan*np.ones(x.shape)  # Initialize with NaN
            value = self._parameter_value_map[var]

            if self._element_ID is None:  # Same constant value for all elements
                variables[var] = value*np.ones(x.shape)
                continue

            # Individual mapping
            indices = np.where(np.isin(self._parameter_value_map['element_ID'], self._element_ID))[0]
            ind_opp = np.where(np.isin(self._element_ID, self._parameter_value_map['element_ID']))[0]
            if len(value)==1:  # Same value for all relevant elements
                variables[var][ind_opp] = value
            else:
                variables[var][ind_opp] = value[indices]

        return variables

