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
# This reader loads an text file with nan-delimited shoreline and/or island polygons
# (typically high-resolution). This is used to perform "on_land" check within the 
# simulation. 
# This reader is expected to be useful for higher-resolution simulation with
# complex shorelines. If high-resolution shoreline is not essential, it is 
# probably more efficient to use reader_global_landmask()
# 
#
# 
# __on_land__() checks are performed with shapely.vectorized.contains()
# 
# see some example below:
# 
# https://jorisvandenbossche.github.io/blog/2017/03/18/vectorized-shapely-cython/
# https://gist.github.com/pelson/9785576
# https://automating-gis-processes.github.io/site/notebooks/L1/geometric-objects.html

from opendrift.readers.basereader import BaseReader, ContinuousReader
import warnings
import pyproj
from shapely.geometry import Polygon, MultiPolygon,asPolygon
import shapely
import shapely.vectorized
import numpy as np

class Reader(BaseReader,ContinuousReader):
    """
    A "custom" landmask reader which can read user-input text file with
    nan-delimited polygons (in geographic WGS84 projection).
    This reader is expected to be useful for higher-resolution simulation with
    complex shorelines. If high-resolution shoreline is not essential, it is
    probably more efficient to use :py:mod:`.reader_global_landmask`.
    .. code::
        NaN NaN
        172.7027815248  -40.5000246431
        172.7041782184  -40.5016429927
        ...
        172.6000000000  -40.5000000000
        172.7027815248  -40.5000246431
        NaN NaN
        175.8076904447  -41.3622783036
        ...
    Args:
        polygon_file: path to file with polyons.
    """

    name = 'custom_landmask'
    return_block = False
    variables = ['land_binary_mask']
    proj4 = None
    crs   = None
    skippoly = False

    def __init__(self,
                 polygon_file = None,
                 proj4_str = '+proj=lonlat +ellps=WGS84'):

        # self.proj4 = '+proj=lonlat +ellps=WGS84'
        self.proj4 = proj4_str 
        self.crs = pyproj.CRS(self.proj4)
        # self.skippoly = skippoly

        super (Reader, self).__init__ ()
        
        shore = np.loadtxt(polygon_file) # nan-delimited closed polygons for land and islands 
        
        if 'lonlat' not in proj4_str: # not functional yet - polygons must be in wgs84 for now   
            print('custom landmask polygon must be in WGS84 coordinates..for now')
            import pdb;pdb.set_trace()
            # if not in wgs84 coordinates already, convert now to lon,lat
            xx,yy = self.xy2lonlat(shore[:,0], shore[:,1])
            xx[np.where(np.isinf(xx))] = np.nan
            yy[np.where(np.isinf(yy))] = np.nan
            # update the shore polygon
            shore[:,0]=xx.copy()
            shore[:,1]=yy.copy()
            # update the projection
            self.proj4 = '+proj=lonlat +ellps=WGS84'
            self.crs = pyproj.CRS(self.proj4)
            # this doesnt work because self.proj.crs.is_geographic remains false in basereader.py 
            # not sure why/how ?
            import pdb;pdb.set_trace()
        
        # Depth
        self.z = None
        self.xmin, self.ymin = np.nanmin(shore[:,0]), np.nanmin(shore[:,1])
        self.xmax, self.ymax = np.nanmax(shore[:,0]), np.nanmax(shore[:,1])
        # convert to xy
        self.xmin, self.ymin = self.lonlat2xy(self.xmin, self.ymin)
        self.xmax, self.ymax = self.lonlat2xy(self.xmax, self.ymax)

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

    def __on_land__(self, x, y):
        """
        Check if particles positions (x,y) are within any shore or island
        polygon.
        """
        return shapely.vectorized.contains(self.mask,x,y)

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

    def plot(self, variable=None, vmin=None, vmax=None,
             filename=None, title=None, buffer=1, lscale='auto'):
        """Plot geographical coverage of reader."""

        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        
        fig = plt.figure()
        plt.ion()

        corners = self.xy2lonlat([self.xmin, self.xmin, self.xmax, self.xmax],
                                 [self.ymax, self.ymin, self.ymax, self.ymin])
        lonmin = np.min(corners[0]) - buffer*2
        lonmax = np.max(corners[0]) + buffer*2
        latmin = np.min(corners[1]) - buffer
        latmax = np.max(corners[1]) + buffer
        latspan = latmax - latmin

        # Initialise map
        if latspan < 90:
            # Stereographic projection centred on domain, if small domain
            x0 = (self.xmin + self.xmax) / 2
            y0 = (self.ymin + self.ymax) / 2
            lon0, lat0 = self.xy2lonlat(x0, y0)
            sp = ccrs.Stereographic(central_longitude=lon0, central_latitude=lat0)
            ax = fig.add_subplot(1, 1, 1, projection=sp)
            corners_stere = sp.transform_points(ccrs.PlateCarree(), np.array(corners[0]), np.array(corners[1]))
        else:
            # Global map if reader domain is large
            sp = ccrs.Mercator()
            ax = fig.add_subplot(1, 1, 1, projection=sp)
            #map = Basemap(np.array(corners[0]).min(), -89,
            #              np.array(corners[0]).max(), 89,
            #              resolution='c', projection='cyl')
        
        if False:
            # dont plot GHSS coastlines here
            # 
            # GSHHS coastlines
            f = cfeature.GSHHSFeature(scale=lscale, levels=[1],
                                      facecolor=cfeature.COLORS['land'])
            ax.add_geometries(
                f.intersecting_geometries([lonmin, lonmax, latmin, latmax]),
                ccrs.PlateCarree(),
                facecolor=cfeature.COLORS['land'],
                edgecolor='black')

        gl = ax.gridlines(ccrs.PlateCarree())
        gl.xlabels_top = False

        # Get boundary
        npoints = 10  # points per side
        x = np.array([])
        y = np.array([])
        x = np.concatenate((x, np.linspace(self.xmin, self.xmax, npoints)))
        y = np.concatenate((y, [self.ymin]*npoints))
        x = np.concatenate((x, [self.xmax]*npoints))
        y = np.concatenate((y, np.linspace(self.ymin, self.ymax, npoints)))
        x = np.concatenate((x, np.linspace(self.xmax, self.xmin, npoints)))
        y = np.concatenate((y, [self.ymax]*npoints))
        x = np.concatenate((x, [self.xmin]*npoints))
        y = np.concatenate((y, np.linspace(self.ymax, self.ymin, npoints)))
        # from x/y vectors create a Patch to be added to map
        lon, lat = self.xy2lonlat(x, y)
        lat[lat>89] = 89.
        lat[lat<-89] = -89.
        p = sp.transform_points(ccrs.PlateCarree(), lon, lat)
        xsp = p[:, 0]
        ysp = p[:, 1]
        
        # add shapely features of custom_landmak
        if variable is None:
            for poly in self.mask.__dict__['context']:
                # to access xy of shapely prepared feature:
                # self.mask.__dict__['context'][0].exterior.xy
                # 
                # The data are defined in lat/lon coordinate system, so PlateCarree() is the appropriate choice
                # https://scitools.org.uk/cartopy/docs/latest/tutorials/understanding_transform.html
                ax.add_geometries([poly], crs=ccrs.PlateCarree(), facecolor='b', edgecolor='red', alpha=0.8)
            # could maybe use feature.ShapelyFeature ?
            # feature.ShapelyFeature(self.mask, crs = ccrs.PlateCarree())

            buf = (xsp.max()-xsp.min())*.1  # Some whitespace around polygon
            buf = 0
            try:
                ax.set_extent([xsp.min()-buf, xsp.max()+buf, ysp.min()-buf, ysp.max()+buf], crs=sp)
            except:
                pass
        if title is None:
            plt.title(self.name)
        else:
            plt.title(title)
        plt.xlabel('Time coverage: %s to %s' %
                   (self.start_time, self.end_time))

        if variable is not None:
            rx = np.array([self.xmin, self.xmax])
            ry = np.array([self.ymin, self.ymax])
            data = self.get_variables_derived(variable, self.start_time,
                                      rx, ry, block=True)
            rx, ry = np.meshgrid(data['x'], data['y'])
            rlon, rlat = self.xy2lonlat(rx, ry)
            #map_x, map_y = map(rlon, rlat, inverse=False)
            data[variable] = np.ma.masked_invalid(data[variable])
            if hasattr(self, 'convolve'):
                from scipy import ndimage
                N = self.convolve
                if isinstance(N, (int, np.integer)):
                    kernel = np.ones((N, N))
                    kernel = kernel/kernel.sum()
                else:
                    kernel = N
                self.logger.debug('Convolving variables with kernel: %s' % kernel)
                data[variable] = ndimage.convolve(
                            data[variable], kernel, mode='nearest')
            #map.pcolormesh(map_x, map_y, data[variable], vmin=vmin, vmax=vmax)
            ax.pcolormesh(rlon, rlat, data[variable], vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
            #cbar = map.colorbar()
            #cbar.set_label(variable)

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except:
            pass

        if filename is not None:
            plt.savefig(filename)
            plt.close()
        else:
            plt.show()