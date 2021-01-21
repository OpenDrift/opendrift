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
# Copyright 2019, Ole Baadshaug, MET Norway & Ron Saper, Carleton University Canada.
#
# Caveat: This copyright will not interfere with the open nature of OpenDrift and OpenBerg

import numpy as np
from datetime import timedelta, datetime
import pyproj

from opendrift.readers.basereader import BaseReader, ContinuousReader

class Reader(BaseReader, ContinuousReader):
	variables = ['x_sea_water_velocity', 'y_sea_water_velocity']

	def __init__ (self, obslon, obslat, obstime, obsfile=None,
					wind_east=0, wind_north=0, windreader=None, wind_factor=0.018,
					name=None):
		"""	Reader which statistically extrapolate current forcing.
	 	It uses the residual track obtained by subtracting the wind forcing component
	 	from the past observed motion of a particle.
		"""
		if name is not None:
			self.name = name
		else:
			self.name = 'reader_current_from_observation'

		# Cover whole earth, no validity radius yet
		self.proj4 = '+proj=latlong'
		self.xmin = -180
		self.xmax = 180
		self.ymin = -90
		self.ymax = 90

		self.z = np.ma.array([0.], mask=[False],fill_value=1e+20)


		self.start_time = obstime[1]+timedelta(hours=1)
		self.end_time = self.start_time+timedelta(days=10)
		self.times = np.arange(self.start_time, self.end_time, timedelta(hours=1)).astype(datetime)
		self.time_step = self.times[-1] - self.times[-2]

		time_delta_seconds = (obstime[1]-obstime[0]).total_seconds()


		if windreader is not None:
			x0,y0 = windreader.lonlat2xy(obslon[0], obslat[0])
			x1,y1 = windreader.lonlat2xy(obslon[1], obslat[1])

			#NB! Use wind speeds at observed positions
			wind0 = windreader.get_variables_interpolated_xy(['x_wind', 'y_wind'], x=x0, y=y0, z = self.z, time=obstime[0])[0]
			wind1 = windreader.get_variables_interpolated_xy(['x_wind', 'y_wind'], x=x1, y=y1, z = self.z, time=obstime[1])[0]
			meanwind_x = (wind0['x_wind'][0]+ wind1['x_wind'][0])/2
			meanwind_y = (wind0['y_wind'][0]+ wind1['y_wind'][0])/2
		else:
			meanwind_x = wind_east
			meanwind_y = wind_north # wind_east/north should be given in m/s


		g = pyproj.Geod(ellps='WGS84')
		self.azimuth, backazimuth, self.dist = \
			g.inv(obslon[0], obslat[0], obslon[1], obslat[1], radians=False)

		self.speed = self.dist/(time_delta_seconds)

		#Calculate the average speed in x/y-direction of the residual
		self.x_sea_water_velocity = np.atleast_1d(self.speed*np.sin(np.radians(self.azimuth)) - meanwind_x*wind_factor)
		self.y_sea_water_velocity = np.atleast_1d(self.speed*np.cos(np.radians(self.azimuth)) - meanwind_y*wind_factor)


		super(Reader, self).__init__()

	def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None):

		requested_variables, time, x, y, z, outside = self.check_arguments(requested_variables, time, x, y, z)
		nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
			self.nearest_time(time)

		variables = {}
		variables['x_sea_water_velocity'] = self.x_sea_water_velocity*np.ones(len(x))
		variables['y_sea_water_velocity'] = self.y_sea_water_velocity*np.ones(len(y))

		return variables

