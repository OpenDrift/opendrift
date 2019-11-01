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
from opendrift_landmask_data import Landmask

import warnings
import pyproj

class Reader(BaseReader):
    """
    The global landmask reader is based on the coastline data from
    GSHHG (https://www.ngdc.noaa.gov/mgg/shorelines/) optimized for
    checking against landmasks.
    """
    name = 'global_landmask'
    return_block = False
    variables = ['land_binary_mask']
    proj4 = None
    crs   = None
    skippoly = False

    def __init__(self,
                 extent = None,
                 llcrnrlon = None,
                 llcrnrlat = None,
                 urcrnrlon = None,
                 urcrnrlat = None,
                 skippoly = False):
        """
        Initialize land mask using GSHHS dataset.

        Args:
            extent (array): minx, miny, maxx, maxy bounding box in source CRS for which to include
                geometries. Default None (all geometries).

            llcrnrlon (float): minx (Deprecated in favor of extent).

            llcrnrlat (float): miny (Deprecated in favor of extent).

            urcrnrlon (float): maxx (Deprecated in favor of extent).

            urcrnrlat (float): maxy (Deprecated in favor of extent).

            skippoly (bool): use only rasterized landmask to determine whether points are on land.
        """

        self.proj4 = '+proj=lonlat +ellps=WGS84'
        self.crs = pyproj.CRS(self.proj4)
        self.skippoly = skippoly

        super (Reader, self).__init__ ()

        # Depth
        self.z = None

        if extent and (llcrnrlat or llcrnrlon or urcrnrlat or urcrnrlon):
            raise Exception ("'extent' cannot be given togheter with any of 'llcrnrlon', ..'")
        elif extent is None and (llcrnrlat or llcrnrlon or urcrnrlat or urcrnrlon):
            warnings.warn("llcrnrlon, llcrnrlat, et. al. is deprecated for the global landmask reader. Prefer 'extent' in stead.", DeprecationWarning)
            extent = [llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat]

        # Read and store min, max and step of x and y
        if extent is not None:
            self.xmin, self.ymin, self.xmax, self.ymax = extent
        else:
            self.xmin, self.ymin = -180, -90
            self.xmax, self.ymax = 180, 90

        self.xmin, self.ymin = self.lonlat2xy(self.xmin, self.ymin)
        self.xmax, self.ymax = self.lonlat2xy(self.xmax, self.ymax)

        # setup landmask
        self.mask = Landmask(extent, skippoly)

    def __on_land__(self, x, y):
        return self.mask.contains (x, y, skippoly = self.skippoly, checkextent = False)

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

