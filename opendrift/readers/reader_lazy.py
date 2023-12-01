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

import logging; logging.captureWarnings(True); logger = logging.getLogger(__name__)
from opendrift.readers.basereader import BaseReader
from opendrift.readers import reader_from_url, reader_netCDF_CF_generic


class Reader:
    '''For lazy initialisation'''

    logger = logging.getLogger('opendrift')  # using common logger

    def __init__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs
        self.initialised = False

        # Interpret if argument is a string (url, filename) or a prepared dataset
        if isinstance(args[0], str):
            self._lazyname = 'LazyReader: ' + args[0]
        else:
            self._lazyname = 'LazyReader: ' + kwargs['name']
            self._dataset = args[0]
        logger.debug('Delaying initialisation of ' + self._lazyname)

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
        logger.debug('Initialising: ' + self._lazyname)
        self.reader = None

        # Zarr is a special case, handle netCDF_CF_generic a prepared dataset and credentials
        if 'zarr_credentials' in self._kwargs:
            logger.debug('Lazy reader seems to be zarr, calling reader_netCDF_CF_generic')
            self.reader = reader_netCDF_CF_generic.Reader(filename=self._dataset, zarr_storage_options=self._kwargs['zarr_credentials'], name=self._lazyname)
        else:
            self.reader = reader_from_url(self._args[0])

        if self.reader is None:
            raise ValueError('Reader could not be initialised')
        else:
            if 'prepare_args' in dir(self):
                self.reader.prepare(**self.prepare_args)
            logger.debug('Reader initialised: ' + self.reader.name)
            self.initialised = True

    def prepare(self, **kwargs):
        self.prepare_args = kwargs

    def __repr__(self):
        if self.initialised is True:
            return self.reader.__repr__()
        else:
            return 'Lazy reader:\n%s\n%s' % (
                        self._args, self._kwargs)
