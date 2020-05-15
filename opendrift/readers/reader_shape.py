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

from opendrift.readers.basereader import BaseReader
from opendrift_landmask_data import Landmask

import pyproj
import shapely
import shapely.ops
import cartopy

class Reader(BaseReader):
    """
    The shape-file reader (.shp) can be used to load generic shape files the 'landmask' variable.

    Args:
        :param shapes: shape-file or files (.shp)
        :type shapes: string or list of strings.

        :param proj4_str: Proj.4 string of shape file projection coordinates
                          (default: '+proj=lonlat +ellps=WGS84').
        :type proj4_str: string.
    """
    name = 'shape'
    return_block = False
    variables = ['land_binary_mask']
    proj4 = None
    crs   = None
    polys = None
    land  = None
    shapes = None
    always_valid = True

    def __init__(self, shapes, proj4_str = '+proj=lonlat +ellps=WGS84'):
        self.proj4 = proj4_str
        self.crs = pyproj.CRS(self.proj4)
        if isinstance(shapes, list):
            self.shapes = shapes
        else:
            self.shapes = [ shapes ]

        super (Reader, self).__init__ ()

        # Depth
        self.z = None

        # reading shapefiles
        self.polys = []
        for shp in self.shapes:
            self.logger.debug("Reading shapefile: %s" % shp)
            from cartopy import io
            reader = io.shapereader.Reader(shp)
            self.polys.extend(reader.geometries())

        assert len(self.polys) > 0, "no geometries loaded"

        self.logger.info("Pre-processing %d geometries" % len(self.polys))
        self.land = shapely.ops.unary_union(self.polys)

        self.xmin, self.ymin, self.xmax, self.ymax = self.land.bounds
        self.xmin, self.ymin = self.lonlat2xy(self.xmin, self.ymin)
        self.xmax, self.ymax = self.lonlat2xy(self.xmax, self.ymax)

    def __on_land__(self, x, y):
        return shapely.vectorized.contains(self.land, x, y)

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
        return { 'x' : x, 'y' : y, 'land_binary_mask': self.__on_land__(x,y) }


