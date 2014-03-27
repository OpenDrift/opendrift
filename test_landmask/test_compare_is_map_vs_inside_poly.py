#!/usr/bin/env python

# Performance test for checking whether points are on land
#
# For two different map areas (regional and global)
# and different GSHHS resolutions (c, i, f)
# it is tested how many seconds it takes to:
# - initialise map
# - check stranding of a large number of points using
#   - Basemap.is_land()
#   - matplotlib.nxutils.points_inside_polys

# See:
# http://stackoverflow.com/questions/13796315/plot-only-on-continent-in-matplotlib

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.nxutils as nx
from datetime import datetime

def points_in_polys(points, polys):
	# Function to quickly check stranding of many points
	insidePoly = np.array([False]*len(points))
	for poly in polys:
		# NB: use contains_points for matplotlib version >= 1.2.0
		insidePoly[nx.points_inside_poly(points, poly)] = True
	return insidePoly

resolutions = ['c', 'i', 'f']

maps = {'global': [-180, 180, -80, 80],
		'regional': [-16, 16, 50, 70]}

maps = {'regional': [4, 7, 59, 61]}
maps = {'regional': [-16, 16, 50, 70]}
resolutions = ['f']

print '########################################'
for mapname, bounds in maps.iteritems():
	xmin, xmax, ymin, ymax = bounds
	print mapname + ' area: ' + str(bounds)
	for resolution in resolutions:
		print '   ' + resolution# + ' coastline resolution'

		# Initialising map
		startTime = datetime.now()
		map = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,
					  urcrnrlon=xmax, urcrnrlat=ymax, 
					  resolution=resolution, projection='cyl')
		map.drawcoastlines()
		map.drawmapboundary(fill_color='aqua')
		map.fillcontinents(color='coral',lake_color='aqua')
		print '\t%6.2f seconds to initialise map' % \
				(datetime.now()-startTime).total_seconds()

		# Extract polygons for faster checking of stranding
		polys = [p.boundary for p in map.landpolygons]
		
		# Check beaching of random points within map bounds
		npoints = 1000
		isLand = [0]*npoints
		np.random.seed(1)
		x = np.random.uniform(xmin, xmax, npoints)
		y = np.random.uniform(ymin, ymax, npoints)

		# using Basemap.is_land()
		startTime = datetime.now()
		isLandBasemap = [map.is_land(X, Y) for X,Y in zip(x,y)]
		print '\t%6.2f seconds to check that %s points out of %s are ' \
				'stranded, using Basemap.is_land()' % \
				((datetime.now()-startTime).total_seconds(), \
				sum(isLandBasemap), npoints)

		# using polygons
		startTime = datetime.now()
		isLandPoly = points_in_polys(np.c_[x, y], polys)
		print '\t%6.2f seconds to check that %s points out of %s are ' \
				'stranded, using matplotlib polygons' % \
				((datetime.now()-startTime).total_seconds(), \
				sum(isLandPoly), npoints)

	print 
print '########################################'

isLandBasemap = np.array(isLandBasemap)

print sum(isLandBasemap)
print sum(isLandPoly)

map.plot(x[isLandPoly], y[isLandPoly], '*y')
map.plot(x[~isLandPoly], y[~isLandPoly], '*k')
plt.savefig('coastline_poly.png')
map.plot(x[isLandBasemap], y[isLandBasemap], '*y')
map.plot(x[~isLandBasemap], y[~isLandBasemap], '*k')
plt.savefig('coastline_is_land.png')
