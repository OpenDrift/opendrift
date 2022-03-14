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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import logging

logger = logging.getLogger(__name__)

# TODO: replace this with weakref + thread-safety
__roaring_mask__ = None

def get_mask_type():
    """
    Returns whether the landmask is provided by the opendrift_landmask_data
    package (0), or the roaring-landmask package (1).
    """
    try:
        from roaring_landmask import RoaringLandmask as _
        return 1
    except ImportError:
        return 0

def get_mask(skippoly = False, extent = None):
    """
    Returns an instance of the landmask type and landmask. The mask data is
    usually shared between threads.
    """
    global __roaring_mask__

    try:
        if __roaring_mask__ is None:
            from roaring_landmask import RoaringLandmask
            __roaring_mask__ = RoaringLandmask.new()

        mask = __roaring_mask__
        mask_type = 1

    except ImportError:
        from opendrift_landmask_data import Landmask
        mask_type = 0
        mask = Landmask(extent, skippoly)

    return mask_type, mask

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

    def show_landmask_old(landmask):
        maxn = 512.
        dx = (lonmax - lonmin) / maxn
        dy = (latmax - latmin) / maxn
        dx = max(landmask.dx, dx)
        dy = max(landmask.dy, dy)

        x = np.array([lonmin, lonmax])
        y = np.array([latmin, latmax])
        ndx = np.int32(dx / landmask.dx)  # normalized delta
        ndy = np.int32(dy / landmask.dy)

        xm, ym = landmask.invtransform * (x, y)
        xm = xm.astype(np.int32)
        ym = ym.astype(np.int32)
        # TODO:
        # this may fail for small maps if xm[0]==xm[1] or ym[0]==ym[1]

        img = landmask.mask[ym[0]:ym[1]:ndy, xm[0]:xm[1]:ndx]

        from matplotlib import colors
        cmap = colors.ListedColormap([ocean_color, land_color])
        ax.imshow(img, origin = 'lower', extent=[lonmin, lonmax, latmin, latmax],
                  zorder=0,
                    transform=ccrs.PlateCarree(globe=globe), cmap=cmap)

    extent = [lonmin, latmin, lonmax, latmax]

    if fast:
        mask_type, mask = get_mask(True, extent)
        if mask_type == 0:
            show_landmask_old(mask)
        else:
            show_landmask_roaring(mask)
    else:
        if get_mask_type() == 0:
            logger.debug("Using existing GSHHS shapes..")
            _, mask = get_mask(False, extent)
            extent = box(lonmin, latmin, lonmax, latmax)
            extent = shapely.prepared.prep(extent)
            polys = [p for p in mask.polys.geoms if extent.intersects(p)]

            ax.add_geometries(polys,
                    ccrs.PlateCarree(globe=globe),
                    zorder=2,
                    facecolor=land_color,
                    edgecolor='black')
        else:
            logger.debug ("Adding GSHHS shapes..")
            f = cfeature.GSHHSFeature(scale=lscale, levels=[1],
                    facecolor=land_color)
            ax.add_geometries(
                    f.intersecting_geometries([lonmin, lonmax, latmin, latmax]),
                    ccrs.PlateCarree(globe=globe),
                    facecolor=land_color,
                    edgecolor='black')



class Reader(BaseReader, ContinuousReader):
    """
    The global landmask reader is based on the coastline data from
    GSHHG (https://www.ngdc.noaa.gov/mgg/shorelines/) optimized for
    checking against landmasks.

    Args:
        :param extent: minx, miny, maxx, maxy bounding box in source CRS for which to include
            geometries. Default None (all geometries).
        :type extent: array of floats, optional

        :param llcrnrlon: minx.
        :type llcrnrlon: float, optional, deprecated in favor of extent.

        :param llcrnrlat: miny.
        :type llcrnrlat: float, optional, deprecated in favor of extent.

        :param urcrnrlon: maxx.
        :type urcrnrlon: float, optional, deprecated in favor of extent.

        :param urcrnrlat: maxy.
        :type urcrnrlat: float, optional, deprecated in favor of extent.

        :param skippoly: use only rasterized landmask to determine whether points are on land.
        :type skippoly: bool, optional, default False.
    """
    name = 'global_landmask'
    variables = ['land_binary_mask']
    proj4 = None
    crs = None
    skippoly = False
    extent = None
    mask_type = 0  # 0 is opendrift_landmask_data, 1 is roaring_landmask. try not to make any code depend on this!

    def __init__(self,
                 extent=None,
                 llcrnrlon=None,
                 llcrnrlat=None,
                 urcrnrlon=None,
                 urcrnrlat=None,
                 skippoly=False):
        self.proj4 = '+proj=lonlat +ellps=WGS84'
        self.crs = pyproj.CRS(self.proj4)
        self.skippoly = skippoly

        super(Reader, self).__init__()

        # Depth
        self.z = None

        if extent and (llcrnrlat or llcrnrlon or urcrnrlat or urcrnrlon):
            raise Exception(
                "'extent' cannot be given togheter with any of 'llcrnrlon', ..'"
            )
        elif extent is None and (llcrnrlat or llcrnrlon or urcrnrlat
                                 or urcrnrlon):
            warnings.warn(
                "llcrnrlon, llcrnrlat, et. al. is deprecated for the global landmask reader. Prefer 'extent' in stead.",
                DeprecationWarning)
            extent = [llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat]

        # setup landmask
        self.mask_type, self.mask = get_mask(skippoly, extent)

        # Read and store min, max and step of x and y
        if self.mask_type == 0 and extent is not None:
            self.xmin, self.ymin, self.xmax, self.ymax = extent
        else:
            self.xmin, self.ymin = -180, -90
            self.xmax, self.ymax = 180, 90

        if self.mask_type == 0 and extent:
            extent = box(*[self.xmin, self.ymin, self.xmax, self.ymax])
            self.extent = shapely.prepared.prep(extent)

    def __on_land__(self, x, y):
        x = self.modulate_longitude(x)

        if self.mask_type == 0:
            return self.mask.contains(x,
                                      y,
                                      skippoly=self.skippoly,
                                      checkextent=False)
        else:
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
