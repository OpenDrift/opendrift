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

from datetime import datetime
import numpy as np

from opendrift.readers.basereader import BaseReader, ContinuousReader


class Reader(BaseReader, ContinuousReader):
    '''Returning values oscillating in time with given amplitude and period'''

    def __init__(self, variable, amplitude, period_seconds=3600*24,
                 phase=0, zero_time=datetime(2017, 1, 1, 0)):
        '''init with '''

        self.variables = [variable]
        self.amplitude = amplitude
        self.period_seconds = period_seconds
        self.zero_time = zero_time
        self.proj4 = '+proj=latlong +datum=WGS84'
        self.xmin = -180
        self.xmax = 180
        self.ymin = -90
        self.ymax = 90
        self.start_time = None
        self.end_time = None
        self.time_step = None
        self.name = 'oscillating_reader'

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, z=None):

        variables = {'time': time, 'x': x, 'y': y, 'z': z}

        phase = ((time - self.zero_time).total_seconds()/
                 self.period_seconds)*np.pi
        value = self.amplitude*np.sin(phase)
        variables[self.variables[0]] = value*np.ones(x.shape)

        return variables
