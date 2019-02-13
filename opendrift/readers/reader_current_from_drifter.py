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
from datetime import datetime, timedelta

import numpy as np
import pyproj

from opendrift.readers.basereader import BaseReader


class Reader(BaseReader):

    return_block = False
    variables = ['x_sea_water_velocity', 'y_sea_water_velocity']

    def __init__(self, lons, lats, times, wind=None, name='reader_current_from_drifter'):
        """Reconstruct currents from time series of drifter positions"""

        self.name = name
        lons = np.asarray(lons)
        lats = np.asarray(lats)

        # Cover whole earth, no validity radius yet
        self.proj4 = '+proj=latlong'
        self.xmin = -180
        self.xmax = 180
        self.ymin = -90
        self.ymax = 90

        self.times = times[0:-1]
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        # NB: assuming her constant time step
        self.time_step = self.times[-1] - self.times[-2]

        lon1 = lons[0:-1]
        lon2 = lons[1::]
        lat1 = lats[0:-1]
        lat2 = lats[1::]
        self.lon = (lon1+lon2)/2
        self.lat = (lat1+lat2)/2

        g = pyproj.Geod(ellps='WGS84')
        self.azimuth, backazimuth, self.dist = \
            g.inv(lon1, lat1, lon2, lat2, radians=False)
        self.speed = []
        print self.dist
        for i, t in enumerate(self.dist):
            self.speed.append(self.dist[i] /
                              (times[i+1] - times[i]).total_seconds())
        self.x_sea_water_velocity = self.speed*np.sin(np.radians(self.azimuth))
        self.y_sea_water_velocity = self.speed*np.cos(np.radians(self.azimuth)) 

        # Subtract wind/stokes contribution
        #self.x_sea_water_velocity = self.x_sea_water_velocity - 0.01*12

        ## Subtract Stokes from Wam 
        #self.stokes_u = np.ones(self.lon.shape)*0
        #self.stokes_v = np.ones(self.lon.shape)*0
        #self.red_fac = np.ones(self.lon.shape)*0
        #from opendrift.readers import reader_grib
        #from opendrift.models.physics_methods import stokes_drift_profile_breivik
        #reader_wam = reader_grib.Reader(
        # '/vol/hindcast3/nofo_olje_paa_vann_2015/WamOILgeo.grb_20150610_0300')

        #for i, t in enumerate(self.times[0:-1]):
        #    if i == 0:
        #        continue
        #    var, prof = reader_wam.get_variables_interpolated(variables=
        #    ['sea_surface_wave_significant_height',
        #     'sea_surface_wave_stokes_drift_y_velocity',
        #     'sea_surface_wave_period_at_variance_spectral_density_maximum',
        #     'sea_surface_wave_stokes_drift_x_velocity'],
        #     time=t, lon=np.atleast_1d(self.lon[i]),
        #             lat=np.atleast_1d(self.lat[i]), block=True, z=np.atleast_1d(0))
        #             #lat=np.atleast_1d(self.lat[i]), block=True, z=np.array(0))
        #    z = -0.64  # Drifter depth
        #    #z = -0.3
        #    wave_peak_period = var[
        #        'sea_surface_wave_period_at_variance_spectral_density_maximum']
        #    self.stokes_u[i], self.stokes_v[i], S = \
        #        stokes_drift_profile_breivik(
        #            var['sea_surface_wave_stokes_drift_x_velocity'],
        #            var['sea_surface_wave_stokes_drift_y_velocity'],
        #            var['sea_surface_wave_significant_height'],
        #            wave_peak_period, z)  # Should be mean period, not peak

        #self.stokes_u[0] = self.stokes_u[1]
        #self.stokes_v[0] = self.stokes_v[1]

        # Run constructor of parent Reader class
        super(Reader, self).__init__()


    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):
        
        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        variables = {}

        try:
            stokes_x = self.stokes_u[indxTime]
            stokes_y = self.stokes_v[indxTime]
        except:
            stokes_x = 0
            stokes_y = 0

        l = np.ones(len(x))
        variables['x_sea_water_velocity'] = \
            self.x_sea_water_velocity[indxTime]*l - stokes_x
        variables['y_sea_water_velocity'] = \
            self.y_sea_water_velocity[indxTime]*l - stokes_y

        return variables
