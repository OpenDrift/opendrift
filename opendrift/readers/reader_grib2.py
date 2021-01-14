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

    def __init__(self, filename, proj4=None, engine='cfgrib'):
        """
        Grib file reader

        Args:

            filename: path to grib or grib2 file.

            proj4: Optional projection override.

            engine: grib engine used by xarray, default is `cfgrib`.

        Returns: Grib-file reader.

        """
        logger.warning("This reader is experimental and may change in breaking ways without a major version change")

        self.dataset = xr.open_dataset(filename, engine=engine)
        self.name = filename

        self.times = np.atleast_1d(self.dataset.time.values)

        # convert to python datetime
        assert '[ns]' in self.times.dtype.str, "datetime is not in nanoseconds"
        self.times = np.array([
            datetime.utcfromtimestamp(d.astype('int') / 1e9)
            for d in self.times
        ])

        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        if len(self.times) > 1:
            self.time_step = self.times[1] - self.times[0]
            assert np.all(
                np.diff(self.times) == self.time_step), "time step not constant"

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
                        assert self.proj4 == v.attrs[
                            'GRIB_gridType'], "all variables must have the same projection"
        else:
            logger.info("Using supplied projection: %s" % self.proj4)

        # parsing variables
        self.variable_mapping = {}

        self.lon = self.dataset.longitude.values
        self.lat = self.dataset.latitude.values

        for v in self.dataset.variables:
            var = self.dataset[v]
            std_name = var.standard_name

            if std_name is not None:
                std_name = self.variable_aliases.get(std_name, std_name)
                self.variable_mapping[std_name] = v
                logger.debug("Found standard variable: %s" % v)

        self.variables = list(self.variable_mapping.keys())

        super().__init__()

        self._make_projected_grid_(self.lon, self.lat)

        self.set_buffer_size(5.)

    def get_variables(self,
                      requested_variables,
                      time=None,
                      x=None,
                      y=None,
                      z=None):

        requested_variables, time, x, y, z, _outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        logger.debug("Request variables: %s" % requested_variables)

        ix0, ix1, iy0, iy1 = self._bbox_(x, y)

        variables = {}

        # This works because we have already asserted that delta_x and delta_y
        # are approx. constant.
        variables['y'] = self.x[ix0:ix1]
        variables['x'] = self.y[iy0:iy1]

        variables['time'] = nearestTime
        variables['z'] = z

        for v in requested_variables:
            par = v
            if hasattr(self, 'rotate_mapping') and par in self.rotate_mapping:
                logger.debug('Using %s to retrieve %s' %
                    (self.rotate_mapping[par], par))
                if par not in self.variable_mapping:
                    self.variable_mapping[par] = \
                        self.variable_mapping[
                            self.rotate_mapping[par]]

            var = self.variable_mapping[v]
            var = self.dataset[var]
            logger.debug("Fetching %s [%d:%d, %d:%d]" %
                         (v, ix0, ix1, iy0, iy1))
            variables[v] = self._slice_variable_(var, indxTime, slice(iy0, iy1), slice(ix0, ix1), 0, 0).values

        return variables

