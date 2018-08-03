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

import logging
from opendrift.readers.basereader import BaseReader
from opendrift.readers import reader_from_url


class Reader(object):
    '''For lazy initialisation'''
    
    def __init__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs
        self.initialised = False
        self._lazyname = 'LazyReader: ' + args[0]
        logging.debug('Delaying initialisation of ' + self._lazyname)

    def __getattr__(self, name):
        if self.initialised is False:
            if name == 'name':
                return self._lazyname
            if name == 'is_lazy':
                return True
            self.initialise()

        try:
            return object.__getattribute__(self.reader, name)
        except:
            return object.__getattribute__(self, name)
        
    def get_variables(self, *args, **kwargs):
        return self.reader.get_variables(*args, **kwargs)

    def initialise(self):
        logging.debug('Initialising: ' + self._lazyname)
        self.reader = reader_from_url(self._args[0])
        if self.reader is None:
            raise ValueError('Reader could not be initialised') 
        else:
            logging.debug('Reader initialised: ' + self.reader.name) 
            self.initialised = True

    def __repr__(self):
        if self.initialised is True:
            return self.reader.__repr__()
        else:
            return 'Lazy reader:\n%s\n%s' % (
                        self._args, self._kwargs)
