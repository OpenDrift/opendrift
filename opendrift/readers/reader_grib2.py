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

import numpy as np
import xarray as xr
from datetime import datetime

from .basereader import BaseReader
from .basereader.structured import StructuredReader

class Reader(BaseReader, StructuredReader):
    dataset = None
    variables = None

    def __init__(self, filename, proj4 = None, engine = 'cfgrib'):
        """
        Grib file reader

        Args:

            filename: path to grib or grib2 file.

            proj4: Optional projection override.

            engine: grib engine used by xarray, default is `cfgrib`.

        Returns: Grib-file reader.

        """
        self.dataset = xr.open_dataset(filename, engine = engine)
        self.name = filename


        self.times = self.dataset.time.values

        # convert to python datetime
        assert '[ns]' in self.times.dtype.str, "datetime is not in nanoseconds"
        self.times = np.array([datetime.utcfromtimestamp(d.astype('int') / 1e9) for d in self.times])

        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        self.time_step = self.times[1] - self.times[0]
        assert np.all(np.diff(self.times) == self.time_step), "time step not constant"

        # scan for projection
        self.proj4 = proj4
        if self.proj4 is None:
            for v in self.dataset.variables:
                v = self.dataset[v]
                if 'GRIB_gridType' in v.attrs:
                    if self.proj4 is None:
                        self.proj4 = v.attrs['GRIB_gridType']
                        logger.info("projection: %s" % self.proj4)
                    else:
                        assert self.proj4 == v.attrs['GRIB_gridType'], "all variables must have the same projection"
        else:
            logger.info("Using supplied projection: %s" % self.proj4)

        # parsing variables
        self.variable_mapping = {}

        for v in self.dataset.variables:
            var = self.dataset[v]
            std_name = var.attrs.get('GRIB_cfName', None)

            # TODO
            if std_name == 'eastward_wind':
                std_name = 'x_wind'

            if std_name is not None:
                self.variable_mapping[std_name] = v
                logger.debug("Found standard variable: %s" % v)

        self.variables = list(self.variable_mapping.keys())

        super().__init__()

        logger.debug("Finding bounds of reader")
        self.X, self.Y = self.lonlat2xy(self.dataset.longitude.values, self.dataset.latitude.values)
        self.xmin, self.xmax = np.min(self.X[:]), np.max(self.X[:])
        self.ymin, self.ymax = np.min(self.Y[:]), np.max(self.Y[:])

        self.delta_x = np.diff(self.X).flat[0]
        self.delta_y = np.diff(self.Y, axis = 0).flat[0]

        eq_eps = 1.e-1

        assert np.all(np.abs(self.delta_x - np.diff(self.X)) < eq_eps), "Grid is not equidistant in X direction"
        assert np.all(np.abs(self.delta_y - np.diff(self.Y, axis = 0)) < eq_eps), "Grid is not equidistant in Y direction"

        self.x = self.X[0,:]
        self.y = self.Y[:,0]

        assert np.all(np.abs(np.tile(self.x, (self.X.shape[0], 1)) - self.X) < eq_eps), "X coordinates are not aligned along Y direction"
        assert np.all(np.abs(np.tile(np.atleast_2d(self.y).T, (1, self.Y.shape[1])) - self.Y) < eq_eps), "Y coordinates are not aligned along X direction"

        self.set_buffer_size(5.)

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None):

        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        if z is not None:
            assert np.all(z == self.dataset.heightAboveGround.item()), "height does not match layer"

        logger.debug("Request variables: %s" % requested_variables)

        ix0, ix1, iy0, iy1 = self._bbox_contains_(x, y, self.buffer)
        print(ix0, ix1, iy0, iy1)

        variables = {}

        # This works because we have already asserted that delta_x and delta_y
        # are approx. constant.
        variables['y'] = self.x[ix0:ix1]
        variables['x'] = self.y[iy0:iy1]

        variables['time'] = nearestTime
        variables['z'] = z

        for v in requested_variables:
            var = self.variable_aliases.get(v, v)
            var = self.variable_mapping[v]
            var = self.dataset[var]
            logger.debug("Fetching %s [%d:%d, %d:%d]" % (v, ix0, ix1, iy0, iy1))
            variables[v] = var[indxTime, 0, ix0:ix1, iy0:iy1].values


        return variables

    def _bbox_contains_(self, x, y, buffer = 0):
        """
        Find bounding box on grid containing points (x, y)

        Args:

            buffer: number of indices to add as buffer
        """
        ix = (x - self.xmin) / self.delta_x
        ix0, ix1 = np.min(ix), np.max(ix)

        iy = (y - self.ymin) / self.delta_y
        iy0, iy1 = np.min(iy), np.max(iy)

        ix0 = np.max((0, ix0 - buffer)).astype(int)
        iy0 = np.max((0, iy0 - buffer)).astype(int)

        ix1 = np.min((len(self.x), ix1 + buffer)).astype(int)
        iy1 = np.min((len(self.y), iy1 + buffer)).astype(int)

        return (ix0, ix1, iy0, iy1)


