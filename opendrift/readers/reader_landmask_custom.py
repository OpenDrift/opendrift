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
# Copyright 2020 Simon Weppe. MetService NZ
# 
# This reader loads an text file with nan-delimited shoreline or island polygons
# (typically high-resolution). This is used to perform "on_land" check within the 
# simulation.
# This reader is expected to be useful for higher-resolution simulation with
# complex shorelines. If high-resolution shoreline is not essential, it is 
# probably more efficient to use reader_global_landmask()
# 
#
# 
# in_polygong checks are performed with shapely.vectorized.contains()
# see some example below:
# 
# https://jorisvandenbossche.github.io/blog/2017/03/18/vectorized-shapely-cython/
# https://gist.github.com/pelson/9785576
# https://automating-gis-processes.github.io/site/notebooks/L1/geometric-objects.html

# check if partciles positions (x,y) are within a shore polygon

from shapely.geometry import Point, Polygon, MultiPolygon,asPolygon
import shapely
import shapely.vectorized
import numpy as np

class Reader(BaseReader):
    """
    The "custom" landmask reader requires a user-input text file with 
    nan-delimited polygons (in geographic WGS84 projection)

    Args:
        :param polygon_file: None.
        :type polygon_file: string.
    """
    name = 'custom_landmask'
    return_block = False
    variables = ['land_binary_mask']
    proj4 = None
    crs   = None
    skippoly = False

    def __init__(self,
                 polygon_file = None):
        self.proj4 = '+proj=lonlat +ellps=WGS84' 
        self.crs = pyproj.CRS(self.proj4)
        # self.skippoly = skippoly

        super (Reader, self).__init__ ()
        
        # shore = np.loadtxt('sydney_shoreline_opendrift.shore') # nan-delimited closed polygons
        shore = np.loadtxt(polygon_file) # nan-delimited closed polygons for land and islands
        # Depth
        self.z = None
        self.xmin, self.ymin = np.min(shore[:,0]), np.min(shore[:,1])
        self.xmax, self.ymax = np.max(shore[:,0]), np.max(shore[:,1])
        # convert to xy ?
        self.xmin, self.ymin = self.lonlat2xy(self.xmin, self.ymin)
        self.xmax, self.ymax = self.lonlat2xy(self.xmax, self.ymax)

        # setup landmask
        # self.mask = Landmask(extent, skippoly)
        
        # Loop through polygon to build MultiPolygon object
        # 
        # make sure that start and end lines are [nan,nan] as well
        if not np.isnan(shore[0,0]):
            shore = np.vstack(([np.nan,np.nan],shore))
        if not np.isnan(shore[-1,0]):
            shore = np.vstack((shore,[np.nan,np.nan]))
        id_nans = np.where(np.isnan(shore[:,0]))[0]
        poly = []
        for cnt,id_i in enumerate(id_nans[:-1]):
            # The shapely.geometry.asShape() family of functions can be used to wrap Numpy coordinate arrays 
            # https://shapely.readthedocs.io/en/latest/manual.html
            poly.append(asPolygon(shore[id_nans[cnt]+1:id_nans[cnt+1]-1,:]))
        # We can pass multiple Polygon -objects into our MultiPolygon as a list
        self.mask= MultiPolygon(poly)
        self.mask = shapely.prepared.prep(self.mask)

        # lon_rel = np.linspace(151.21,151.215,100) # start from second grid point rather than actual edge
        # lat_rel = np.linspace(-33.86,-33.856,100)
        # [lon_grid,lat_grid]=np.meshgrid(lon_rel,lat_rel)
        # x_tmp = np.ravel(lon_grid).copy()
        # y_tmp = np.ravel(lat_grid).copy()

        # on_land = shapely.vectorized.contains(land, np.ravel(lon_grid), np.ravel(lat_grid))

        # import matplotlib.pyplot as plt 
        # plt.ion()
        # # plot mesh nodes
        # plt.plot(x_tmp,y_tmp,'b.')
        # plt.plot(x_tmp[on_land],y_tmp[on_land],'r.')
        # plt.plot(shore[:,0],shore[:,1],'k')
        # import pdb;pdb.set_trace()

    def __on_land__(self, x, y):
        # return self.mask.contains (x, y, skippoly = self.skippoly, checkextent = False)
        return shapely.vectorized.contains(self.mask,x,y)
        # return a boolean array , True if in land/island poly, False if not

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

        >> check format is correct 1/0 and not boolean ??