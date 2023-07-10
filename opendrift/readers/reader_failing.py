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
# Copyright 2023, Knut-Frode Dagestad, MET Norway

from opendrift.readers.basereader import BaseReader, ContinuousReader


class Reader(BaseReader, ContinuousReader):
    '''A reader for testing discarding after failure'''

    def __init__(self):

        self.variables = ['x_wind', 'y_wind']
        self.proj4 = '+proj=latlong'
        self.xmin = -180
        self.xmax = 180
        self.ymin = -90
        self.ymax = 90
        self.start_time = None
        self.end_time = None
        self.time_step = None
        self.name = 'failing_reader'

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, z=None):

        raise ValueError('Failing reader, for testing only.')

