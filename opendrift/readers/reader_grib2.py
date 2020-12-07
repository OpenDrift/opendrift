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
# Copyright 2020, Gaute Hope, MET Norway

import logging
logger = logging.getLogger(__name__)

import xarray as xr

from .basereader import BaseReader
from .basereader.structured import StructuredReader

class Reader(BaseReader, StructuredReader):
    dataset = None

    def __init__(self, filename, engine = 'cfgrib'):
        """
        Grib file reader

        Args:

            filename: path to grib or grib2 file.

            engine: grib engine used by xarray, default is `cfgrib`.

        Returns: Grib-file reader.

        """
        self.dataset = xr.open_dataset(filename, engine = engine)

        self.xmin = self.dataset.longitude.min()
        self.xmax = self.dataset.longitude.max()
        self.ymin = self.dataset.latitude.min()
        self.ymax = self.dataset.latitude.max()

        super().__init__()

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None):
        pass


