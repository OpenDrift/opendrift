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
    def from_shpfiles(shpfiles, proj4_str = '+proj=lonlat +ellps=WGS84', invert=False):
        """
        Construct a shape-reader from shape-files (.shp)

        Args:
            :param shapes: shape-file or files (.shp)
            :type shapes: string or list of file names as strings.

            :param proj4_str: Proj.4 string of shape file projection coordinates
                            (default: '+proj=lonlat +ellps=WGS84').
            :type proj4_str: string.
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
            shp_iters.append(reader.records())
        return Reader(itertools.chain(*shp_iters), proj4_str, invert=invert)

    def _ingest_landmask_only(self):
        assert len(self.polys) > 0, "no geometries loaded"
        logger.info("Pre-processing %d geometries" % len(self.polys))
        self.land = {
            self.variables[0]: shapely.ops.unary_union(self.polys)
        }

    def _ingest_multivariable(self):
        self.land = {}
        for var in self.variables:
            features = filter(lambda x: x.attributes['variable'] == var,
                self.records)
            self.land[var] = shapely.ops.unary_union([x.geometry
                for x in features])

    def __init__(self, shapes, proj4_str = '+proj=lonlat +ellps=WGS84', invert=False):
        self.invert = invert  # True if polygons are lakes and not land areas
        self.proj4 = proj4_str
        self.crs = pyproj.CRS(self.proj4)
        super (Reader, self).__init__ ()
        # Depth
        self.z = None
        try:
            # Test if shapes are records or geometries
            self.records = list(shapes)
            self.polys = [x.geometry for x in self.records]
        except AttributeError:
            # This is a classic land mask file
            self.polys = list(shapes) # reads geometries if Shape Readers
            self.records = None
            self._ingest_landmask_only()
        else:
            try:
                # When records are passed, check for the `variable` attribute
                add_vars = set([x.attributes['variable'] for x in self.records])
            except KeyError:
                # It's a plain land mask file without attributes
                self._ingest_landmask_only()
            else:
                # It's a multi-variable shapefile
                self.variables = list(add_vars)
                self._ingest_multivariable()
        # Regardless of the number of variables, the extent is always the sum
        # of the polygons
        whole = shapely.ops.unary_union(self.polys)
        self.xmin, self.ymin, self.xmax, self.ymax = whole.bounds
        self.xmin, self.ymin = self.lonlat2xy(self.xmin, self.ymin)
        self.xmax, self.ymax = self.lonlat2xy(self.xmax, self.ymax)

    def __on_land__(self, x, y, var):
        if self.invert is False:
            return shapely.vectorized.contains(self.land[var], x, y)
        else:  # Inverse if polygons are lakes and not land areas
            return 1 - shapely.vectorized.contains(self.land[var], x, y)

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
        requestedVariables = \
            self.check_arguments(requestedVariables, time, x, y, z)[0]
        variables = { 'x' : x, 'y' : y}
        for var in requestedVariables:
            variables[var] = self.__on_land__(x, y, var)
        return variables

