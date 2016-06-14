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

import logging
import numpy as np
from readers.reader import pyproj


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
        try:
            self.config['drift']['scheme']
        except:
            self.config['drift']['scheme'] = 'euler'
        # Runge-Kutta scheme
        if self.config['drift']['scheme'] == 'runge-kutta':
            x_vel = self.environment.x_sea_water_velocity
            y_vel = self.environment.y_sea_water_velocity
            # Calculate x,y from lon,lat
            start_x, start_y = self.lonlat2xy(self.elements.lon,
                                              self.elements.lat)
            # Find midpoint
            az = np.degrees(np.arctan2(x_vel, y_vel))
            speed = np.sqrt(x_vel*x_vel + y_vel*y_vel)
            dist = speed*self.time_step.total_seconds()*.5
            geod = pyproj.Geod(ellps='WGS84')
            mid_lon, mid_lat, dummy = geod.fwd(self.elements.lon,
                                               self.elements.lat,
                                               az, dist, radians=False)
            # Find current at midpoint, a half timestep later
            logging.debug('Runge-kutta, fetching half time-step later...')
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

    def advect_wind(self, wind_drift_factor=None):
        # Elements at ocean surface (z=0) are advected with given percentage
        # of wind speed. NB: Only basic Euler schema is implemented

        if wind_drift_factor == None:
            wind_drift_factor = self.elements.wind_drift_factor

        # Convert wind_drift_factor to array
        if len(np.atleast_1d(wind_drift_factor)) == 1:
            wind_drift_factor = wind_drift_factor*np.ones(len(self.elements))
        # Convert z to array
        if not hasattr(self.elements, 'z'):
            z = np.zeros(len(self.elements))  # Assumed on surface, if no z
        else:
            z = self.elements.z
        if len(np.atleast_1d(z)) == 1:
            z = z*np.ones(len(self.elements))

        wind_drift_factor[self.elements.z < 0] = 0

        x_wind = self.environment.x_wind
        y_wind = self.environment.y_wind
        try:
            if self.config['drift']['relative_wind'] == True:
                # Use wind relative to (subtracted) ocean current
                x_wind = x_wind - self.environment.x_sea_water_velocity
                y_wind = y_wind - self.environment.y_sea_water_velocity
        except:
            pass

        self.update_positions(x_wind*wind_drift_factor,
                              y_wind*wind_drift_factor)

    def stokes_drift(self):

        if np.max(np.array(
            self.environment.sea_surface_wave_stokes_drift_x_velocity)) == 0:
            logging.debug('No Stokes drift velocity available.')
            return

        logging.debug('Calculating Stokes drift')
        # Assuming deep water:
        wave_peak_period = self.environment.\
            sea_surface_wave_period_at_variance_spectral_density_maximum
        wavelength = 9.81/(2*3.14)*np.power(wave_peak_period, 2)
        depth_reduction_factor = np.exp(4*3.14*self.elements.z/wavelength)

        self.update_positions(
            self.environment.sea_surface_wave_stokes_drift_x_velocity*
                depth_reduction_factor,
            self.environment.sea_surface_wave_stokes_drift_y_velocity*
                depth_reduction_factor)

    def wave_mixing(self, time_step_seconds):
        return  # To be implemented by subclasses, e.g. downward mixing of oil

    def resurface_elements(self, minimum_depth):
        # Keep surfacing elements in water column as default,
        # i.e. no formation of surface slick
        surface = np.where(self.elements.z >= 0)[0]
        self.elements.z[surface] = minimum_depth  

    def deactivate_stranded_elements(self):
        self.deactivate_elements(self.environment.land_binary_mask == 1,
                                 reason='stranded')

    def wind_speed(self):
        return np.sqrt(self.environment.x_wind**2 +
                       self.environment.y_wind**2)

    def significant_wave_height(self):
        # Significant wave height, parameterise from wind if not available
        if hasattr(self.environment, 'sea_surface_wave_significant_height'):
            return self.environment.sea_surface_wave_significant_height
        else:
            # Neumann and Pierson, 1966
            return 0.2*np.power(self.wind_speed,2)/9.81

    def wave_frequency(self):
        return 2*np.pi/self.environment. \
            sea_surface_wave_period_at_variance_spectral_density_maximum

    def wave_period(self):
        return self.environment.\
                sea_surface_wave_period_at_variance_spectral_density_maximum

    def wave_energy(self):
        return 9.81*1028*np.power(self.significant_wave_height(),2)/16
    
    def wave_energy_dissipation(self):
        # Delvigne and Sweeney
        return 0.0034*self.sea_water_density()*9.81*\
            np.power(self.significant_wave_height(), 2)

    def wave_damping_coefficient(self):
        return (10E-5)*2*np.pi*self.wave_frequency()* \
            np.power(self.wave_energy(), 0.25)

    #def sea_water_density(self):
    #    return 1027  # kg/m3

    def sea_surface_wave_breaking_fraction(self):
        f = 0.032*(self.wind_speed() - 5)/self.wave_period()
        f[f<0] = 0
        return f

    def air_density(self):
        return 1.225  # Could calculate from temperature

    def windspeed_from_stress(self):
        wind_stress = np.sqrt(self.environment.surface_downward_x_stress**2+
                              self.environment.surface_downward_y_stress**2)
        return windspeed_from_stress_polyfit(wind_stress)

def wind_drag_coefficient(windspeed):
    '''Large and Pond (1981), J. Phys. Oceanog., 11, 324-336.'''
    Cd = 0.0012*np.ones(len(windspeed))
    Cd[windspeed>11] = 0.001*(0.49 + 0.065*windspeed[windspeed>11])
    return Cd

def windspeed_from_stress_polyfit(wind_stress):
    '''Inverting Large and Pond (1981) using polyfit'''
    windspeed = np.linspace(0, 30, 30)
    rho_air = 1.225
    stress = wind_drag_coefficient(windspeed)*rho_air*(windspeed**2)
    z = np.polyfit(stress, windspeed, 3)
    p = np.poly1d(z)
    return p(wind_stress)


