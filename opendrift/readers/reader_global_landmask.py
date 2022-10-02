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

import warnings
import pyproj
import numpy as np
import shapely.vectorized
import shapely.prepared
from shapely.geometry import box
import shapely.wkb as wkb
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import logging

logger = logging.getLogger(__name__)

# TODO: replace this with weakref + thread-safety
__roaring_mask__ = None
__polys__ = None

def get_mask():
    """
    Returns an instance of the landmask type and landmask. The mask data is
    usually shared between threads.
    """
    global __roaring_mask__

    if __roaring_mask__ is None:
        from roaring_landmask import RoaringLandmask
        __roaring_mask__ = RoaringLandmask.new()

    return __roaring_mask__

class LandmaskFeature(cfeature.GSHHSFeature):
    def __init__(self, scale='auto', globe=None, **kwargs):
        super().__init__(scale, **kwargs)

        if globe is not None:
            self._crs = ccrs.PlateCarree(globe=globe)

    def geometries(self):
        self.intersecting_geometries(extent=None)

    def intersecting_geometries(self, extent):
        global __polys__

        if self._scale == 'auto':
            scale = self._scale_from_extent(extent)
        else:
            scale = self._scale[0]

        # # If scale is full use the geometries from roaring landmask, otherwise
        # # fall back to Cartopy provider.
        # if scale == 'f':
        #     logger.debug(f"Adding full GSHHG shapes from roaring-landmask, extent: {extent}..")
        #     if __polys__ is None:
        #         from roaring_landmask import Gshhg
        #         polys = Gshhg.wkb()
        #         __polys__ = wkb.loads(polys)

        #     if extent is not None:
        #         extent = box(extent[0], extent[2], extent[1], extent[3])
        #         extent = shapely.prepared.prep(extent)
        #         return (p for p in __polys__.geoms if extent.intersects(p))
        #     else:
        #         return __polys__.geoms
        # else:
        logger.debug(f"Adding GSHHG shapes from cartopy, scale: {scale}, extent: {extent}..")
        return super().intersecting_geometries(extent)

def plot_land(ax, lonmin, latmin, lonmax, latmax, fast, ocean_color = 'white', land_color = cfeature.COLORS['land'], lscale = 'auto', globe=None):
    """
    Plot the landmask or the shapes from GSHHG.
    """
    def show_landmask_roaring(roaring):
        maxn = 512.
        dx = (lonmax - lonmin) / maxn
        dy = (latmax - latmin) / maxn
        dx = max(roaring.dx, dx)
        dy = max(roaring.dy, dy)

        x = np.arange(lonmin, lonmax, dx)
        y = np.arange(latmin, latmax, dy)

        yy, xx = np.meshgrid(y, x)
        img = roaring.mask.contains_many_par(xx.ravel(), yy.ravel()).reshape(yy.shape).T

        from matplotlib import colors
        cmap = colors.ListedColormap([ocean_color, land_color])
        ax.imshow(img, origin = 'lower', extent=[lonmin, lonmax, latmin, latmax],
                  zorder=0,
                    transform=ccrs.PlateCarree(globe=globe), cmap=cmap)
    if fast:
        show_landmask_roaring(get_mask())
    else:
        land = LandmaskFeature(scale=lscale, facecolor=land_color, globe=globe)

        ax.add_feature(land, zorder=2,
                       facecolor=land_color,
                       edgecolor='black')



class Reader(BaseReader, ContinuousReader):
    """
    The global landmask reader is based on the coastline data from
    GSHHG (https://www.ngdc.noaa.gov/mgg/shorelines/) optimized for
    checking against landmasks.
    """
    name = 'global_landmask'
    variables = ['land_binary_mask']
    proj4 = None
    crs = None

    def __init__(self):
        self.proj4 = '+proj=lonlat +ellps=WGS84'
        self.crs = pyproj.CRS(self.proj4)

        super(Reader, self).__init__()

        # Depth
        self.z = None

        self.xmin, self.ymin = -180, -90
        self.xmax, self.ymax = 180, 90

        # setup landmask
        self.mask = get_mask()

    def __on_land__(self, x, y):
        x = self.modulate_longitude(x)
        x = x.astype(np.float64)
        y = y.astype(np.float64)
        return self.mask.contains_many_par(x, y)

    def get_variables(self,
                      requestedVariables,
                      time=None,
                      x=None,
                      y=None,
                      z=None):
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
        return {'land_binary_mask': self.__on_land__(x, y)}
