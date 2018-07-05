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
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

from opendrift.readers.basereader import BaseReader


class Reader(object):
    '''For lazy initialisation'''
    
    def __init__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs
        self.initialised = False
        print('Delaying initialisation')

    def __getattr__(self, name):
        if self.initialised is False:
            print('Initialising...')
            from opendrift.models.oceandrift import OceanDrift
            o = OceanDrift()
            self.reader = o._reader_fom_url(self._args[0])
            if self.reader is None:
                raise ValueError('Reader could not be initialised') 
            else:
                self.initialised = True
        return object.__getattribute__(self.reader, name)
        
    def get_variables(self, *args, **kwargs):
        pass

    def __repr__(self):
        return 'Lazy reader:\n%s\n%s' % (
                     self._args, self._kwargs)
