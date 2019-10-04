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
# Copyright 2019, Gaute Hope, MET Norway

from opendrift.readers.basereader import BaseReader
import opendrift_landmask_data as old

import logging
import warnings
import shapely.wkb
import shapely.vectorized
import cartopy
import numpy as np
import os.path

class Reader(BaseReader):
    name = 'cartopy_landmask'
    return_block = False
    variables = ['land_binary_mask']
    proj4 = None

    __land_shapes__ = None
    __land__ = None
    _crs_    = None

    def __init__(self, resolution = 'c'):
        """
        Initialize cartopy land mask using GSHHS dataset.

        Args:

            resolution (str): Resolution category for landmask. Default is
                'c'. Possible values: 'c', 'l', 'i', 'h', 'f'.
        """
        reader = cartopy.feature.GSHHSFeature(scale = resolution, level = 2) # only used for CRS

        # The cartopy reader CRS is a unit-radian scaled earth, so that lon, lat is equivalent, but
        # projection conversion between full scale ellipsoids do not work. For that it is better to use
        # PlateCarree with a WGS84 globe.
        self._crs_ = reader.crs
        self.proj4 = self._crs_.proj4_init

        super (Reader, self).__init__ ()

        # Depth
        self.z = None

        # Time
        self.start_time = None
        self.end_time = None
        self.time_step = None

        # Read and store min, max and step of x and y
        self.xmin, self.ymin = -180, -90
        self.xmax, self.ymax = 180, 90

        self.xmin, self.ymin = self.lonlat2xy(self.xmin, self.ymin)
        self.xmax, self.ymax = self.lonlat2xy(self.xmax, self.ymax)
        self.delta_x = None
        self.delta_y = None

        resmap = { 'c': 0, 'l': 1, 'i': 2, 'h': 3, 'f': 4 }
        cachef = old.GSHHS[resmap[resolution]]

        logging.debug("cartopy: loading cache from: %s", cachef)
        with open(cachef, 'rb') as fd:
            self.__land_shapes__ = shapely.wkb.loads(fd.read())

        self.__land__ = shapely.prepared.prep(self.__land_shapes__)

    def zoom_map(self, buffer=0.2,
                 lonmin=None, lonmax=None, latmin=None, latmax=None):
        logging.warning('Zooming not implemented for cartopy reader: (%s to %s E), (%s to %s N)' %
                     (lonmin, lonmax, latmin, latmax) )

    def __on_land__(self, x, y):
        assert isinstance(x, np.ndarray) and isinstance(y, np.ndarray)

        return np.array(shapely.vectorized.contains(self.__land__, x, y),
                        dtype=np.bool)

    def get_variables(self, requestedVariables, time = None,
                      x = None, y = None, z = None, block = False):
        """
        Get binary mask of whether elements are on land or not.

        Args:
            x (deg[]): longitude (decimal degrees)
            y (deg[]): latitude (decimal degrees)
            ...

        x, y is given in reader local projection.

        Returns:
            Binary mask of point x, y on land.

        """

        self.check_arguments(requestedVariables, time, x, y, z)
        return { 'land_binary_mask': self.__on_land__(x,y) }

