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

from opendrift.readers.basereader import BaseReader, ContinuousReader

import pyproj
import shapely
import shapely.ops
import shapely.vectorized
import cartopy
import itertools
import logging
import numpy as np
from scipy.spatial import cKDTree
import shapely as shp
from typing import Iterable

logger = logging.getLogger(__name__)

class Reader(BaseReader, ContinuousReader):
    """
    The shape reader can be used to load generic shapes as the 'landmask' variable.

    Args:
        :param shapes: shapely geometries.
        :type shapes: iterable.

        :param proj4_str: Proj.4 string of shape file projection coordinates
                          (default: '+proj=lonlat +ellps=WGS84').
        :type proj4_str: string.
    """
    name = 'shape'
    variables = ['land_binary_mask']
    proj4 = None
    crs   = None
    polys = None
    land  = None
    always_valid = True

    @staticmethod
    def from_shpfiles(shpfiles, proj4_str = '+proj=lonlat +ellps=WGS84', invert=False, land_buffer_distance=0):
        """
        Construct a shape-reader from shape-files (.shp)

        Args:
            :param shapes: shape-file or files (.shp)
            :type shapes: string or list of file names as strings.

            :param proj4_str: Proj.4 string of shape file projection coordinates
                            (default: '+proj=lonlat +ellps=WGS84').
            :type proj4_str: string.
            :type land_buffer_distance: float. Buffer distance around polygons
        """
        if isinstance(shpfiles, list):
            shpfiles = shpfiles
        else:
            shpfiles = [ shpfiles ]

        # reading shapefiles
        shp_iters = []
        for shp in shpfiles:
            logger.debug("Reading shapefile: %s" % shp)
            from cartopy import io
            reader = io.shapereader.Reader(shp)
            shp_iters.append(reader.geometries())

        return Reader(itertools.chain(*shp_iters), proj4_str, invert=invert, land_buffer_distance=land_buffer_distance)

    def __init__(self, shapes, proj4_str = '+proj=lonlat +ellps=WGS84', invert=False, land_buffer_distance=None):

        self._land_kdtree = None
        self._land_kdtree_buffer_distance = land_buffer_distance
        self.invert = invert  # True if polygons are lakes and not land areas
        self.proj4 = proj4_str
        self.crs = pyproj.CRS(self.proj4)

        super (Reader, self).__init__ ()

        # Depth
        self.z = None

        # loading and pre-processing shapes
        self.polys = list(shapes) # reads geometries if Shape Readers
        assert len(self.polys) > 0, "no geometries loaded"

        logger.info("Pre-processing %d geometries" % len(self.polys))
        self.land = shapely.ops.unary_union(self.polys)

        self.xmin, self.ymin, self.xmax, self.ymax = self.land.bounds
        self.xmin, self.ymin = self.lonlat2xy(self.xmin, self.ymin)
        self.xmax, self.ymax = self.lonlat2xy(self.xmax, self.ymax)

    def __on_land__(self, x, y):
        if self.invert is False:
            return shapely.vectorized.contains(self.land, x, y)
        else:  # Inverse if polygons are lakes and not land areas
            return 1 - shapely.vectorized.contains(self.land, x, y)

    def get_variables(self, requestedVariables, time = None,
                      x = None, y = None, z = None):
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
        return { 'x' : x, 'y' : y, 'land_binary_mask': self.__on_land__(x,y) }

    def get_nearest_outside(self, x, y, buffer_distance: float):
        """
        Determine the nearest point outside the loaded polygons, including
        an additional buffer distance
            Args:
                x (deg[]): longitude (decimal degrees) or projected x coordinate
                y (deg[]): latitude (decimal degrees) or projected y coordinate
                buffer_distance: buffer zone around polygons
                ...

            x, y, buffer are given in reader local projection.

            Returns:
                lon, lat, and distance to nearest points

        """
        if self._land_kdtree is None or self._land_kdtree_buffer_distance != buffer_distance:
            logger.info("Building KDTree from %d geometries with buffer distance %e" % (len(self.polys), buffer_distance))

            _land_polygon_buffered = shp.unary_union(self.polys).buffer(buffer_distance)
            xy = unwrap(_land_polygon_buffered)

            self._land_kdtree = cKDTree(xy.transpose())
            self._land_kdtree_buffer_distance = buffer_distance

        dist, ind = self._land_kdtree.query(np.stack((x, y), axis=1))

        return self._land_kdtree.data[ind, 0], self._land_kdtree.data[ind, 1], dist


def unwrap(geom):
    """
    Unwrap a shapely geometry or a list thereof into boundary coordinates

    Returns:
        array of shape (2, N)
    """
    if isinstance(geom, shp.LineString):
        return np.array(geom.xy)
    elif isinstance(geom, shp.MultiPolygon):
        return np.concatenate([unwrap(p.boundary) for p in geom.geoms], axis=1)
    elif isinstance(geom, shp.Polygon):
        return unwrap(geom.boundary)
    elif isinstance(geom, shp.MultiLineString):
        return np.concatenate([unwrap(p) for p in geom.geoms], axis=1)
    elif isinstance(geom, Iterable):
        return np.concatenate([unwrap(p) for p in geom], axis=1)
    else:
        raise ValueError("Unknown type %d" % type(geom))
