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
# Copyright 2015,2020 Knut-Frode Dagestad, MET Norway

from opendrift.readers.basereader import BaseReader, StructuredReader


class Reader(BaseReader, StructuredReader):
    '''Reader based on static 2D arrays of variables'''

    def __init__(self, x, y, array_dict, proj4='+proj=latlong'):
        '''init with a map {'variable_name': array, ...}'''

        self.proj4 = proj4
        self.xmin = x.min()
        self.xmax = x.max()
        self.ymin = y.min()
        self.ymax = y.max()

        self.variables = list(array_dict.keys())
        self.array_dict = {'x': x, 'y': y, 'z': 0}
        for var in array_dict:
            self.array_dict[var] = array_dict[var]

        self.start_time = None
        self.end_time = None
        self.time_step = None
        self.name = 'reader_constant_2d'

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, z=None):

        self.array_dict['time'] = time
        return self.array_dict
