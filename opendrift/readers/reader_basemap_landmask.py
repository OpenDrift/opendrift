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
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import logging

import os
import gc
import numpy as np
import collections
import matplotlib
if ('DISPLAY' not in os.environ and
	    'PYCHARM_HOSTED' not in os.environ and
        os.name != 'nt'):
    logging.info('No display found. Using non-interactive Agg backend')
    matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
try:
    import matplotlib.nxutils as nx
    has_nxutils = True  # For matplotlib version < 1.2
except:
    from matplotlib.path import Path
    has_nxutils = False  # For matplotlib version >= 1.2
from matplotlib.collections import PolyCollection
from matplotlib.backends.backend_agg import FigureCanvasAgg

from opendrift.readers.basereader import BaseReader

RasterizedBasemap = collections.namedtuple('RasterizedBasemap', 'xmin xmax ymin ymax resolution data')

class Reader(BaseReader):


    name = 'basemap_landmask'
    return_block = False  # Vector based, so checks only individual points

    # Variables (CF standard names) which
    # can be provided by this model/reader
    variables = ['land_binary_mask']

    def __init__(self, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
                 resolution='i', projection='merc',
                 rasterize=True, rasterize_resolution=500):

        logging.debug('Creating Basemap...')

        if projection == 'cyl':  # TODO: should be fixed
            raise ValueError('Cyl projection not supported, as rasterisation needs units of m')

        # Generate Basemap instane
        self.map = Basemap(llcrnrlon, llcrnrlat,
                           urcrnrlon, urcrnrlat, area_thresh=0,
                           resolution=resolution, projection=projection)

        # Store proj4 string
        self.proj4 = self.map.proj4string

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

        # Depth
        self.z = None

        # Time
        self.start_time = None
        self.end_time = None
        self.time_step = None

        # Read and store min, max and step of x and y
        self.xmin, self.ymin = self.lonlat2xy(llcrnrlon, llcrnrlat)
        self.xmax, self.ymax = self.lonlat2xy(urcrnrlon, urcrnrlat)
        self.delta_x = None
        self.delta_y = None

        # Extract polygons for faster checking of stranding
        if has_nxutils is True:
            self.polygons = [p.boundary for p in self.map.landpolygons]
        else:
            self.polygons = [Path(p.boundary) for p in self.map.landpolygons]
            
        # Generate rasterized version of polygons for faster checking of stranding
        # (if enabled)
        if (rasterize == True):
            if (len(self.map.landpolygons) == 0):
                logging.debug('No land polygons to rasterize...')
                self.bmap_raster = None
            else:
                logging.debug('Creating rasterized Basemap...')
                try:
                    self.bmap_raster = self.gen_land_bitmap(self.map, rasterize_resolution)
                except Exception as e:
                    logging.warning('Rasterizing Basemap failed! Continuing without rasterized version. Received "' + e.message + '"')
                    self.bmap_raster = None
        else:
            self.bmap_raster = None

        # Calculate aspect ratio, to minimise whitespace on figures
        # Drawback is that empty figure is created in interactive mode
        meanlat = (llcrnrlat + urcrnrlat)/2
        aspect_ratio = np.float(urcrnrlat - llcrnrlat) / \
            (np.float(urcrnrlon-llcrnrlon))
        if projection != 'cyl':
            aspect_ratio = aspect_ratio / np.cos(np.radians(meanlat))
        if aspect_ratio > 1:
            self.figsize=(10./aspect_ratio, 10.)
            plt.figure(0, figsize=(10./aspect_ratio, 10.))
        else:
            self.figsize=(11., 11.*aspect_ratio)
            plt.figure(0, figsize=(11., 11.*aspect_ratio))
        ax = plt.axes([.05, .05, .85, .9])
            
    def on_land_polycheck(self, x, y):
        points = np.c_[x, y]
        land = np.zeros_like(x, dtype=np.bool)
        if has_nxutils is True:
            for polygon in self.polygons:
                land[nx.points_inside_poly(points, polygon)] = True
        else:
            for polygon in self.polygons:
                land = land + polygon.contains_points(points)
        return land

    """
    Returns a vector of booleans with True if the point (x[i], y[i]) is on land
    """
    def on_land(self, x, y):
        #return [self.map.is_land(x0, y0) for x0,y0 in zip(x,y)]  # uncomment for simulation in lakes
        x0 = (x - self.bmap_raster.xmin)/self.bmap_raster.resolution
        y0 = (y - self.bmap_raster.ymin)/self.bmap_raster.resolution
        x0 = x0.astype(np.int32)
        y0 = y0.astype(np.int32)
        y0 = self.bmap_raster.data.shape[0] - y0 - 1
        
        #Clip out of bounds
        np.clip(x0, 0, self.bmap_raster.data.shape[1]-1, out=x0)
        np.clip(y0, 0, self.bmap_raster.data.shape[0]-1, out=y0)
            
        land = (self.bmap_raster.data[y0, x0] == 0)
        coords = np.flatnonzero(land)
        logging.debug('Checking ' + str(len(coords)) + ' of ' + str(len(x)) + ' coordinates to polygons')
        
        if (len(coords) > 0):
            land[coords] = self.on_land_polycheck(x[coords], y[coords])
        
        return land

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, z=None, block=False):

        if isinstance(requestedVariables, str):
            requestedVariables = [requestedVariables]

        self.check_arguments(requestedVariables, time, x, y, z)

        # Apparently it is necessary to first convert from x,y to lon,lat
        # using proj library, and back to x,y using Basemap instance
        # Perhaps a bug in Basemap related to x_0/y_0 offsets?
        # Nevertheless, seems not to affect performance
        lon, lat = self.xy2lonlat(x, y)
        x, y = self.map(lon, lat, inverse=False)
        
        if (self.bmap_raster == None):
            insidePoly = self.on_land_polycheck(x, y)
        else:
            insidePoly = self.on_land(x, y)

        variables = {}
        variables['land_binary_mask'] = insidePoly

        return variables

    """
    Rasterizes a Basemap object into a bitmap with a given resolution 
    (each cell has a size of resolution_meters x resolution_meters)
    """
    @staticmethod
    def gen_land_bitmap(bmap, resolution_meters):
                
        #Get land polygons and bbox of polygons
        polys = []
        xmin = np.finfo(np.float64).max
        xmax = -np.finfo(np.float64).max
        ymin = xmin
        ymax = xmax
                
        logging.debug('Rasterizing Basemap, number of land polys: ' + str(len(bmap.landpolygons)))
        # If no polys: return a zero map
        if (len(bmap.landpolygons) == 0):
            raise Exception('Basemap contains no land polys to rasterize')
        
        for polygon in bmap.landpolygons:
            coords = polygon.get_coords()
            xmin = min(xmin, np.min(coords[:,0]))
            xmax = max(xmax, np.max(coords[:,0]))
            ymin = min(ymin, np.min(coords[:,1]))
            ymax = max(ymax, np.max(coords[:,1]))
            polys.append(coords)
            
        xmin = np.floor(xmin/resolution_meters)*resolution_meters
        xmax = np.ceil(xmax/resolution_meters)*resolution_meters
        ymin = np.floor(ymin/resolution_meters)*resolution_meters
        ymax = np.ceil(ymax/resolution_meters)*resolution_meters
        
        # For debugging
        logging.debug('Rasterizing Basemap, bounding box: ' + str([xmin, xmax, ymin, ymax]))
        
        # Switch backend to prevent creating an empty figure in notebook
        orig_backend = plt.get_backend()
        plt.switch_backend('agg')
        
        # Create figure to help rasterize
        fig = plt.figure(frameon=False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)   
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        
        # Set aspect and resolution
        # Aspect gives 1 in high plot
        aspect = (xmax-xmin)/(ymax-ymin)
        resolution_dpi = (ymax-ymin) / resolution_meters
        while resolution_dpi > 10000:
            logging.debug('Too large dpi %s, reducing by factor of 2'
                          % resolution_dpi)
            resolution_dpi = resolution_dpi / 2
        
        fig.set_dpi(resolution_dpi)
        fig.set_size_inches(aspect, 1)
        
        # Add polygons
        lc = PolyCollection(polys, facecolor='k', lw=0)
        ax.add_collection(lc)
        
        # Create canvas and rasterize
        canvas = FigureCanvasAgg(fig)
        try:
            canvas.draw()
            width, height = canvas.get_width_height()
            rgb_data = np.frombuffer(canvas.tostring_rgb(), dtype='uint8').reshape(height, width, 3)
            data = rgb_data[:,:,1]
            plt.close(fig) #comment this for debugging purposes and replace with plt.show()
            logging.debug('Rasterized size: ' + str([width, height]))
        except MemoryError:
            gc.collect()
            raise Exception('Basemap rasterized size too large: ' 
                            + str(aspect*resolution_dpi) + '*' + str(resolution_dpi) 
                            + ' cells')
        finally:
            # Reset backend
            plt.switch_backend(orig_backend)
        
        
        return RasterizedBasemap(xmin, xmax, ymin, ymax, resolution_meters, data)
