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
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader
import shapely.geometry as sgeom
from shapely import wkb, clip_by_rect
import logging
import roaring_landmask

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

        if self._scale == 'auto':
            scale = self._scale_from_extent(extent)
        else:
            scale = self._scale[0]

        if extent is not None:
            extent_geom = sgeom.box(extent[0], extent[2],
                                    extent[1], extent[3])
        for level in self._levels:
            geoms = cfeature.GSHHSFeature._geometries_cache.get((scale, level, extent_geom))
            if geoms is None:
                if scale in ['f', 'full']: # Getting fullres polygons from roaring_landmask
                    geoms = cfeature.GSHHSFeature._geometries_cache.get('ROARING')
                    if geoms is None:
                        logger.debug('Getting fullres shapes from roaring landmask')
                        provider = roaring_landmask.LandmaskProvider.Gshhg
                        shapes = roaring_landmask.Shapes.new(provider)
                        polys = roaring_landmask.Shapes.wkb(provider)
                        __polys__ = wkb.loads(polys)
                        geoms = __polys__.geoms
                        cfeature.GSHHSFeature._geometries_cache['ROARING'] = geoms
                    else:
                        logger.debug('Getting fullres shapes from roaring landmask cache')
                else:
                    logger.debug(f'Loading shapes (\'{scale}\' level {level}) with Cartopy shapereader...')
                    # Load GSHHS geometries from appropriate shape file.
                    # TODO selective load based on bbox of each geom in file.
                    path = shapereader.gshhs(scale, level)
                    geoms = tuple(shapereader.Reader(path).geometries())

                # Unlike Caropy.GSHHSFeature, we clip the polygons to plot extent
                #geoms = [g.intersection(extent_geom) for g in geoms if g.intersects(extent_geom)]
                delta = .1  # Adding small delta to avoid segments along map borders
                geoms = [clip_by_rect(g, extent[0]-delta, extent[2]-delta, extent[1]+delta, extent[3]+delta)
                         for g in geoms if g.intersects(extent_geom)]

                # Due to the above, the cache is also depending on extent
                # Note: if plotting several extents within same session, it would be benefitial
                #       with another cache independent of extent
                cfeature.GSHHSFeature._geometries_cache[(scale, level, extent_geom)] = geoms

            for geom in geoms:
                yield geom

            if scale in ['f', 'full']:
                return  # Only a single (combined) level with full resolution troaring_landmask

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

def plot_land(ax, lonmin, latmin, lonmax, latmax, fast, ocean_color = 'white', land_color = cfeature.COLORS['land'], lscale = 'auto', crs_plot=None, crs_lonlat=None):
    """
    Plot the landmask or the shapes from GSHHG.
    """
    def show_landmask_roaring(roaring):
        maxn = 512.

        if crs_plot is not None and crs_lonlat is not None:
            # Make transformer and convert lon,lat to Mercator x,y
            transformer = pyproj.Transformer.from_crs(crs_lonlat, crs_plot)
            xmin, ymin = transformer.transform(lonmin, latmin)
            xmax, ymax = transformer.transform(lonmax, latmax)
            dx = (xmax-xmin) / maxn
            dy = (ymax-ymin) / maxn
            dx = max(roaring.dx, dx)
            dy = max(roaring.dy, dy)

            x = np.arange(xmin, xmax, dx)
            y = np.arange(ymin, ymax, dy)
            yy, xx = np.meshgrid(y, x)
            lons, lats = transformer.transform(xx, yy, direction='inverse')

            extent=[xmin, xmax, ymin, ymax]
            transform = crs_plot
        else:
            dlon = (lonmax - lonmin) / maxn
            dlat = (latmax - latmin) / maxn
            dlon = max(roaring.dx, dlon)
            dlat = max(roaring.dy, dlat)

            lons = np.arange(lonmin, lonmax, dlon)
            lats = np.arange(latmin, latmax, dlat)
            lats, lons = np.meshgrid(lats, lons)

            extent=[lonmin, lonmax, latmin, latmax]
            transform = crs_lonlat

        img = roaring.mask.contains_many(lons.ravel(), lats.ravel()).reshape(lons.shape).T

        from matplotlib import colors
        cmap = colors.ListedColormap([ocean_color, land_color])

        ax.imshow(img, origin = 'lower',
                  extent=extent, zorder=0, cmap=cmap,
                  transform=transform)

    if fast:
        show_landmask_roaring(get_mask())
    else:
        if latmax-latmin < 3 and lonmax-lonmin < 3:
            # Check first if there is land at all within map bounds
            roaring = get_mask()
            maxn = 512.
            dlon = (lonmax - lonmin) / maxn
            dlat = (latmax - latmin) / maxn
            dlon = max(roaring.dx, dlon)
            dlat = max(roaring.dy, dlat)
            lons = np.arange(lonmin, lonmax, dlon)
            lats = np.arange(latmin, latmax, dlat)
            lats, lons = np.meshgrid(lats, lons)
            land = roaring.mask.contains_many(lons.ravel(), lats.ravel())
            if not land.any():
                logger.debug('No raster land within plot area, skipping plotting of GSHHG vectors')
                return

        land = LandmaskFeature(scale=lscale, facecolor=land_color, globe=crs_lonlat.globe, levels=[1,5,6])

        ax.add_feature(land, zorder=0,
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
        return self.mask.contains_many(x, y)

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
