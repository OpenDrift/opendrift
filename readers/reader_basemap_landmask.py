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

from datetime import datetime, timedelta
import logging

import numpy as np
from mpl_toolkits.basemap import Basemap
try:
    import matplotlib.nxutils as nx
    has_nxutils = True  # For matplotlib version < 1.2
except:
    from matplotlib.path import Path
    has_nxutils = False  # For matplotlib version >= 1.2

from readers import Reader


class Reader(Reader):

    name = 'basemap_landmask'
    return_block = False  # Vector based, so checks only individual points

    # Variables (CF standard names) which
    # can be provided by this model/reader
    variables = ['land_binary_mask']

    def __init__(self, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
                 resolution='i', projection='cyl'):

        logging.debug('Creating Basemap...')

        # Set up Basemap
        self.map = Basemap(llcrnrlon, llcrnrlat,
                           urcrnrlon, urcrnrlat,
                           resolution=resolution, projection=projection)

        # Store proj4 string
        self.proj4 = self.map.proj4string

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

        # Depth
        self.depths = None

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

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, depth=None, block=False):

        if isinstance(requestedVariables, str):
            requestedVariables = [requestedVariables]

        self.check_arguments(requestedVariables, time, x, y, depth)

        # Apparently it is necessary to first convert from x,y to lon,lat
        # using proj library, and back to x,y using Basemap instance
        # Perhaps a bug in Basemap related to x_0/y_0 offsets?
        # Nevertheless, seems not to affect performance
        lon, lat = self.xy2lonlat(x, y)
        x, y = self.map(lon, lat, inverse=False)
        points = np.c_[x, y]

        insidePoly = np.array([False]*len(x))
        if has_nxutils is True:
            for polygon in self.polygons:
                insidePoly[nx.points_inside_poly(points, polygon)] = True
        else:
            for polygon in self.polygons:
                insidePoly += np.array(polygon.contains_points(points))

        variables = {}
        variables['land_binary_mask'] = insidePoly

        return variables
