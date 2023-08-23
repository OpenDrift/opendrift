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
# Copyright 2015, Knut-Frode Dagestad, MET Norway

from datetime import datetime, timedelta

import numpy as np
import pyproj

from opendrift.readers.basereader import BaseReader, ContinuousReader


class Reader(BaseReader, ContinuousReader):
    variables = ['x_sea_water_velocity', 'y_sea_water_velocity']

    def __init__(self, lons, lats, times, wind=None, waves=None, z=0,
                 name='reader_current_from_drifter'):
        """Reconstruct currents from time series of drifter positions"""

        self.name = name
        lons = np.asarray(lons)
        lats = np.asarray(lats)
        valid = np.isfinite(lons+lats)
        if np.sum(valid) < len(lons):
            print('Skipping %s invalid positions' % (len(lons)-np.sum(valid)))
            lons = lons[valid]
            lats = lats[valid]
            times = times[valid]

        # Cover whole earth, no validity radius yet
        self.proj4 = '+proj=latlong'
        self.xmin = -180
        self.xmax = 180
        self.ymin = -90
        self.ymax = 90

        self.times = times[0:-1]
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        # NB: assuming here constant time step
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
        for i, t in enumerate(self.dist):
            self.speed.append(self.dist[i] /
                              (times[i+1] - times[i]).total_seconds())
        self.speed = np.array(self.speed)
        self.x_sea_water_velocity = self.speed*np.sin(np.radians(self.azimuth))
        self.y_sea_water_velocity = self.speed*np.cos(np.radians(self.azimuth))

        # Subtract wind/stokes contribution
        if waves is not None:
            ts = waves.get_timeseries_at_position(lon=self.lon, lat=self.lat,
                    #start_time=self.start_time, end_time=self.end_time,
                    times=self.times,
                    variables=['sea_surface_wave_significant_height',
                               'sea_surface_wave_stokes_drift_x_velocity',
                               'sea_surface_wave_stokes_drift_y_velocity',
                               'sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment'])
            # NB: Stokes drift already rotated to lat-lon axes

            if z <= 0:
                print('Calculating Stokes drift at depth of %s m' % -z)
                from opendrift.models.physics_methods import stokes_drift_profile_monochromatic

                wave_peak_period = ts[
                    'sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment']
                    #'sea_surface_wave_period_at_variance_spectral_density_maximum']
                self.stokes_u, self.stokes_v, S = stokes_drift_profile_monochromatic(
                    ts['sea_surface_wave_stokes_drift_x_velocity'],
                    ts['sea_surface_wave_stokes_drift_y_velocity'],
                    ts['sea_surface_wave_significant_height'],
                    wave_peak_period, z)  # Should be mean period, not peak

                print('Subtracting Stokes drift from current')
                stokes_speed = np.sqrt(self.stokes_u*self.stokes_u+self.stokes_v*self.stokes_v)
                print('Current in range %s - %s' % (self.speed.min(), self.speed.max()))
                print('Stokesdrift in range %s - %s' % (stokes_speed.min(), stokes_speed.max()))
                print('Stokesdrift fraction in range %s - %s' % (
                        (stokes_speed/self.speed).min(), (stokes_speed/self.speed).max()))
                self.x_sea_water_velocity = self.x_sea_water_velocity - self.stokes_u
                self.y_sea_water_velocity = self.y_sea_water_velocity - self.stokes_v

        # Run constructor of parent Reader class
        super(Reader, self).__init__()


    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None):

        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        variables = {}

        l = np.ones(len(x))
        variables['x_sea_water_velocity'] = \
            self.x_sea_water_velocity[indxTime]*l
        variables['y_sea_water_velocity'] = \
            self.y_sea_water_velocity[indxTime]*l

        return variables
