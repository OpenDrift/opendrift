#!/usr/bin/env python

import sys
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import rcParams
from datetime import datetime


def propagate_towards_land(proj, res, xmin, xmax, ymin, ymax, lon, lat):
		rcParams['figure.subplot.hspace'] = 0.1 # less height between subplots
		fig = plt.figure()
		# Initialise map
		startTime = datetime.now()
		map = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,
					  urcrnrlon=xmax, urcrnrlat=ymax, 
					  resolution=res, projection=proj)
		ax = fig.add_subplot(211)
		map.drawcoastlines()
		map.drawmapboundary(fill_color='aqua')
		map.fillcontinents(color='coral',lake_color='aqua')
		initTime = datetime.now()-startTime

		# Propagate single particle until hitting coast
		lonHistory = np.array([lon])
		latHistory = np.array([lat])
		startTime = datetime.now()
		while not map.is_land(lon, lat):
			lon = lon + 0.01
			lat = lat + 0.001
			lonHistory = np.append(lonHistory, lon)
			latHistory = np.append(latHistory, lat)
		calcTime = datetime.now()-startTime

		# Plot and save figure
		plt.text(xmin, ymax, 'initialisation time: ' + 
						str(initTime.total_seconds()) + 's\n'
					'calculation time: ' + 
						str(calcTime.total_seconds()) + 's\n'
					'iterations: ' + str(len(latHistory)) +
					',  resolution: ' + res)
		map.plot(lonHistory, latHistory, '.b')
		map.plot(lon, lat, '*y')
		filename = 'figures_test_landmask_map_size/%s-%s-%s-%s-%s-%s.png' \
				% (proj, res, xmin, xmax, ymin, ymax)

		# Fullres zoom around landing point
		ax = fig.add_subplot(212)
		subsetmap = Basemap(llcrnrlon=lon-0.5,llcrnrlat=lat-.5,
					  urcrnrlon=lon+0.5, urcrnrlat=lat+.5, 
					  resolution='f', projection=proj)
		subsetmap.drawcoastlines()
		subsetmap.drawmapboundary(fill_color='aqua')
		subsetmap.fillcontinents(color='coral',lake_color='aqua')
		subsetmap.plot(lonHistory, latHistory, '*b')
		subsetmap.plot(lon, lat, '.y')
		
		#plt.show()
		plt.savefig(filename)
		print 'Saved figure: ' + filename

# Run tests

# starting point
lon = -2
lat = 59

# Medium sized area, North Europe
propagate_towards_land('cyl', 'c', -16.,16.,50.,70., lon, lat)
propagate_towards_land('cyl', 'i', -16.,16.,50.,70., lon, lat)
propagate_towards_land('cyl', 'h', -16.,16.,50.,70., lon, lat)

# (near) Global area
propagate_towards_land('cyl', 'c', -180.,180.,-80.,80., lon, lat)
propagate_towards_land('cyl', 'i', -180.,180.,-80.,80., lon, lat)
propagate_towards_land('cyl', 'h', -180.,180.,-80.,80., lon, lat)
