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

import logging

from opendrift.readers.basereader import BaseReader

import shapely.vectorized
from shapely.ops import unary_union
import cartopy
import numpy as np


class Reader(BaseReader):
    name = 'cartopy_landmask'
    return_block = False
    variables = ['land_binary_mask']
    proj4: str = None

    __prep__: bool = None
    __land__ = None
    _crs_: cartopy.crs.CRS = None

    def __init__(self,
                 prepare: bool = False,
                 source: str = 'gshhs',
                 resolution: str = None,
                 extent = None,
                 llcrnrlon = None,
                 llcrnrlat = None,
                 urcrnrlon = None,
                 urcrnrlat = None):
        """
        Initialize cartopy land mask

        Args:

            prepare (bool): Optimize landmask by joining all geometries. Speeds up
                later checks, but slows down initialization (default: False).

            source (str): Source for landmask. Either 'gshhs' or 'naturalearth'.
                Default is 'gshhs'.

            resolution (str): Resolution category for given source. Default is
                'c' for gshhs and '10m' for 'naturalearth'. For 'gshhs' this is
                'c', 'l', 'i', 'h', 'f'. For 'naturalearth' this is '110m', '50m',
                '10m'. GSHHS generally provides higher resolution.

            extent (array): minx, miny, maxx, maxy bounding box in source CRS for which to include
                geometries. Default None (all geometries).

            llcrnrlon (float): minx (Deprecated in favor of extent).

            llcrnrlat (float): miny (Deprecated in favor of extent).

            urcrnrlon (float): maxx (Deprecated in favor of extent).

            urcrnrlat (float): maxy (Deprecated in favor of extent).

        """
        print ("setting up cartopy")
        self.__prep__ = prepare
        if source == 'naturalearth':
            if resolution is None:
                resolution = '10m'

            reader = cartopy.feature.NaturalEarthFeature('physical', 'coastline', resolution)

        elif source == 'gshhs':
            if resolution is None:
                resolution = 'c'
            reader = cartopy.feature.GSHHSFeature(scale = resolution, level = 2)

        else:
            raise Exception ("invalid source for cartopy landmask")

        # The cartopy reader CRS is a unit-radian scaled earth, so that lon, lat is equivalent, but
        # projection conversion between full scale ellipsoids do not work. For that it is better to use
        # PlateCarree with a WGS84 globe.
        self._crs_ = reader.crs
        self.proj4 = self._crs_.proj4_init

        if extent and (llcrnrlat or llcrnrlon or urcrnrlat or urcrnrlon):
            raise Exception ("'extent' cannot be given togheter with any of 'llcrnrlon', ..'")
        elif extent is None and (llcrnrlat or llcrnrlon or urcrnrlat or urcrnrlon):
            import warnings
            warnings.warn("llcrnrlon, llcrnrlat, et. al. is deprecated for the cartopy reader. Prefer 'extent' in stead.", DeprecationWarning)
            extent = [llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat]

        super (Reader, self).__init__ ()

        # Depth
        self.z = None

        # Time
        self.start_time = None
        self.end_time = None
        self.time_step = None

        # Read and store min, max and step of x and y
        if extent is not None:
            self.xmin, self.ymin, self.xmax, self.ymax = extent
        else:
            self.xmin, self.ymin = -180, -90
            self.xmax, self.ymax = 180, 90

        self.xmin, self.ymin = self.lonlat2xy(self.xmin, self.ymin)
        self.xmax, self.ymax = self.lonlat2xy(self.xmax, self.ymax)
        print ("cartopy setup, xmin, xmax:", self.xmin, self.xmax)
        self.delta_x = None
        self.delta_y = None

        if self.__prep__:
            # XXX: Can be further optimized by saving a cached version
            land_geom = unary_union(list(reader.intersecting_geometries(extent)))
            self.__land__ = shapely.prepared.prep(land_geom)
        else:
            self.__land__ = list(reader.intersecting_geometries(extent))
        print ("setup cartopy")

    def __on_land__(self, x, y):
        assert isinstance(x, np.ndarray) and isinstance(y, np.ndarray)

        if self.__prep__:
            return np.array(shapely.vectorized.contains(self.__land__, x, y),
                            dtype=np.bool)
        else:
            S = x.shape
            x = x.ravel()
            y = y.ravel()

            r = np.zeros(len(x), dtype=np.bool)
            for l in self.__land__:
                if np.all(r): break
                w = np.logical_not(r)
                r[w] = shapely.vectorized.contains(l, x[w], y[w])

            return r.reshape(S)

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

        print ("cart: x,y=", x,y)
        self.check_arguments(requestedVariables, time, x, y, z)
        return { 'land_binary_mask': self.__on_land__(x,y) }

