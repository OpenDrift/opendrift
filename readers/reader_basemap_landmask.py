#!/usr/bin/env python

from datetime import datetime, timedelta

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.nxutils as nx

from readers import Reader


def points_in_polys(points, polys):
    # Function to quickly check stranding of many points
    insidePoly = np.array([False]*len(points))
    for poly in polys:
        # NB: use contains_points for matplotlib version >= 1.2.0
        insidePoly[nx.points_inside_poly(points, poly)] = True
    return insidePoly


class Reader(Reader):

    # Variables (CF standard names) which
    # can be provided by this model/reader
    variables = ['land_binary_mask']

    def __init__(self, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
                 resolution='i', projection='cyl'):

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
        self.startTime = None
        self.endTime = None
        self.timeStep = None

        # Read and store min, max and step of x and y
        self.xmin, self.ymin = self.lonlat2xy(llcrnrlon, llcrnrlat)
        self.xmax, self.ymax = self.lonlat2xy(urcrnrlon, urcrnrlat)
        self.delta_x = None
        self.delta_y = None

        # Extract polygons for faster checking of stranding
        self.polys = [p.boundary for p in self.map.landpolygons]

    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, depth=None):

        if isinstance(requestedVariables, str):
            requestedVariables = [requestedVariables]

        self.check_arguments(requestedVariables, time, x, y, depth)

        x, y = self.xy2lonlat(x, y)
        isLand = points_in_polys(np.c_[x, y], self.polys)

        return isLand
