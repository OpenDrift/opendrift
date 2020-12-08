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

        self.lat = self.dataset.latitude.values
        self.lon = self.dataset.longitude.values
        self.delta_x = np.abs(self.lon[0,1] - self.lon[0,0]).item()
        print(self.delta_x)

        self.xmin = self.dataset.longitude.min().item()
        self.xmax = self.dataset.longitude.max().item()
        self.ymin = self.dataset.latitude.min().item()
        self.ymax = self.dataset.latitude.max().item()

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
            if self.proj4 == 'fakeproj':
                self.proj4 = None

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

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None):

        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        if z is not None:
            assert np.all(z == self.dataset.heightAboveGround.item()), "height does not match layer"

        logger.debug("Request variables: %s" % requested_variables)
        print("x=", x)
        print("y=", y)

        variables = {}

        ix0, ix1, iy0, iy1 = self._bbox_contains_(self.lon, self.lat, x, y, 1.)

        variables['y'] = self.lat[ix0:ix1, iy0:iy1]
        variables['x'] = self.lon[ix0:ix1, iy0:iy1]
        variables['time'] = nearestTime
        variables['z'] = z

        print(variables)

        for v in requested_variables:
            var = self.variable_aliases.get(v, v)
            var = self.variable_mapping[v]
            var = self.dataset[var]
            logger.debug("Fetching %s [%d:%d, %d:%d]" % (v, ix0, ix1, iy0, iy1))
            variables[v] = var[indxTime, 0, ix0:ix1, iy0:iy1].values


        return variables

    def _bbox_contains_(self, X, Y, x, y, margin = 0):
        """
        Find a bounding box in 2D arrays `X` (e.g. `longitudes`) and `Y` (e.g.
        `latitudes`) where all of `x` and `y` are contained.

        The 2D arrays of `X` and `Y` tend to be ireguarily distanced in
        lon-lat space, but approximately so in the (unknown) projection on a
        surface.

        Args:
            X: 2D array of x-coordinate at grid points

            Y: 2D array of y-coordinate at grid points

            x: 1D array of x-coordinate that must be contained in bbox.

            y: 1D array of y-coordinate that must be contained in bbox.

            margin: extend bbox with margin.

        Returns:
            (x0, x1, y0, y1)

            Start and end indices in 2D arrays `X` and `Y` of rectangular
            bounding box.

        If no box is found, an assertion error is raised.
        """
        xmin = x.min() - margin
        xmax = x.max() + margin
        ymin = y.min() - margin
        ymax = y.max() + margin

        W = np.argwhere(
                np.logical_and(
                    np.logical_and(Y>=ymin, Y<=ymax),
                    np.logical_and(X>=xmin, X<=xmax)))

        assert len(W) > 0, "all points are outside X and Y"

        ix0 = np.min(W[:,0])
        ix1 = np.max(W[:,0])
        iy0 = np.min(W[:,1])
        iy1 = np.max(W[:,1])

        return (ix0, ix1, iy0, iy1)



