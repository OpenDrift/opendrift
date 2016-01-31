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
# Copyright 2016, Knut-Frode Dagestad, MET Norway

import numpy as np


class PhysicsMethods(object):
    """Physics methods to be inherited by OpenDriftSimulation class"""

    @staticmethod
    def sea_water_density(T=10., S=35.):
        '''The function gives the density of seawater at one atmosphere
        pressure as given in :

        N.P. Fofonoff and R.C. Millard Jr.,1983,
        Unesco technical papers in marine science no. 44.

        S   = Salinity in promille of the seawater
        T   = Temperature of the seawater in degrees Celsius
        '''

        R4 = 4.8314E-04
        DR350 = 28.106331

        # Pure water density at atmospheric pressure
        # Bigg P.H. (1967) BR. J. Applied Physics pp.:521-537

        R1 = ((((6.536332E-09 * T - 1.120083E-06) * T + 1.001685E-04) *
              T - 9.095290E-03) * T + 6.793952E-02) * T - 28.263737

        #Seawater density at atmospheric pressure
        #coefficients involving salinity :

        R2 = (((5.3875E-09 * T - 8.2467E-07) * T + 7.6438E-05) *
              T - 4.0899E-03) * T + 8.24493E-01

        R3 = (-1.6546E-06*T+1.0227E-04)*T-5.72466E-03

        #International one-atmosphere equation of state of seawater :

        SIG = R1 + (R4*S + R3*np.sqrt(S) + R2)*S
        Dens0 = SIG + DR350 + 1000.
        return Dens0

    def advect_ocean_current(self):
        # Runge-Kutta scheme
        if self.config['drift']['scheme'] == 'runge-kutta':
            x_vel = self.environment.x_sea_water_velocity
            y_vel = self.environment.y_sea_water_velocity
            # Calculate x,y from lon,lat
            start_x, start_y = self.lonlat2xy(self.elements.lon,
                                              self.elements.lat)
            # Find midpoint
            mid_x = start_x + x_vel*self.time_step.total_seconds()*.5
            mid_y = start_y + y_vel*self.time_step.total_seconds()*.5
            mid_lon, mid_lat = self.xy2lonlat(mid_x, mid_y)
            # Find current at midpoint, a half timestep later
            mid_env, profiles, missing = self.get_environment(
                ['x_sea_water_velocity', 'y_sea_water_velocity'],
                self.time + self.time_step/2,
                mid_lon, mid_lat, self.elements.z, profiles=None)
            # Move particles using runge-kutta velocity
            self.update_positions(mid_env['x_sea_water_velocity'],
                                  mid_env['y_sea_water_velocity'])
        elif self.config['drift']['scheme'] == 'euler':
            # Euler scheme
            self.update_positions(self.environment.x_sea_water_velocity,
                                  self.environment.y_sea_water_velocity)
        else:
            raise ValueError('Drift scheme not recognised: ' +
                              self.config['drift']['scheme'])
