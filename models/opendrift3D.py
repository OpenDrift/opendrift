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

import os
import numpy as np
import logging
from datetime import datetime

from opendrift import OpenDriftSimulation
from elements import LagrangianArray


# Defining the oil element properties
class Lagrangian3DArray(LagrangianArray):
    """Extending LagrangianArray to a passive particle that moves in 3dimensions
    The Particle may be buoyant and/or subject to vertical mixing
    buoyant bahaviour is described by terminal velocity
    """

    variables = LagrangianArray.add_variables([
        ('terminal_velocity', {'dtype': np.float32,
                     'units': 'm/s',
                     'default': 0.})])

class OpenDrift3DSimulation(OpenDriftSimulation):
    """Open source buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway
        
        Generic module for particles that move in 3 dimensions 
        and may be to vertical turbulent mixing
        with the possibility for positive or negative buoyancy

        Particles could be pollution (e.g. oil droplets), plankton, nutrients or sediments

        Under construction.
    """

    ElementType = Lagrangian3DArray

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'sea_surface_wave_significant_height',
                          'sea_surface_wave_to_direction',
                          'sea_ice_area_fraction',
                          'x_wind', 'y_wind', 'land_binary_mask',
                          'sea_floor_depth_below_sea_level',
                          'ocean_vertical_diffusivity',
                          'sea_water_temperature', 'sea_water_salinity', 
                          'sea_water_molecular_viscosity',
                          'upward_sea_water_velocity'
                         ]

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       'sea_surface_wave_to_direction': 0,
                       'sea_ice_area_fraction': 0,
                       'x_wind': 0, 'y_wind': 0,
                       'sea_floor_depth_below_sea_level': 100,
                       'ocean_vertical_diffusivity': 0.02, # m2 s-1
                       'sea_water_temperature': 10.,    
                       'sea_water_salinity': 34.,
                       'sea_water_molecular_viscosity': 0.0014, # kg m-1 s-1
                       'upward_sea_water_velocity': 0.
                      }

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'missing_data': 'gray', 'stranded': 'red',
                     'evaporated': 'yellow', 'dispersed': 'magenta'}

    # Configuration
    configspec = '''
        [input]
            readers = list(min=1, default=list(''))
            [[spill]]
                longitude = float(min=-360, max=360, default=5)
                latitude = float(min=-90, max=90, default=60)
                time = string(default='%s')
        [processes]
            turbulentmixing = boolean(default=True)
            verticaladvection = boolean(default=False)
            dispersion = boolean(default=False)
            diffusion = boolean(default=False)
            evaporation = boolean(default=False)
            emulsification = boolean(default=False)
        [drift]
            wind_drift_factor = float(min=0, max=1, default=0.00)
            current_uncertainty = float(min=0, max=5, default=.1)
            wind_uncertainty = float(min=0, max=5, default=1)
            relative_wind = boolean(default=True)
    ''' % (datetime.now().strftime('%Y-%d-%m %H:00'))

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(OpenDrift3DSimulation, self).__init__(*args, **kwargs)

    def update_terminal_velocity(self):
        """Calculate terminal velocity due to bouyancy from own properties
        and environmental variables
        overload this function to create particle-specific behaviour
        """
        self.elements.terminal_velocity = 0.005
        
    def update(self):
        """Update positions and properties of buoyant particles."""

        self.elements.age_seconds += self.time_step.total_seconds()

        # Calculate windspeed, which is needed for emulsification,
        # and dispersion (if wave height is not given)
        # wind speed can also be used to parameterize eddy diffusivity
        windspeed = np.sqrt(self.environment.x_wind**2 +
                            self.environment.y_wind**2)

        ###################
        # Turbulent Mixing
        ###################
        if self.config['processes']['turbulentmixing'] is True:

            self.update_terminal_velocity() 
            self.vertical_mixing()

        ###################
        # Vertical advection
        ###################
        if self.config['processes']['verticaladvection'] is True:

            self.vertical_advection(self.environment.upward_sea_water_velocity)

        ####################################
        # Deactivate elements hitting land
        ####################################
        self.deactivate_elements(self.environment.land_binary_mask == 1,
                                 reason='stranded')

        ##############################################
        # Simply move particles with ambient current
        ##############################################
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        ##############
        # Wind drag
        ##############
        # wind drag should only be activated for particles that are very near the surface
        wind_drift_factor = self.config['drift']['wind_drift_factor']
        if self.config['drift']['relative_wind'] is True:
            self.update_positions((self.environment.x_wind -
                                   self.environment.x_sea_water_velocity)*
                                   wind_drift_factor,
                                  (self.environment.y_wind -
                                   self.environment.y_sea_water_velocity)*
                                   wind_drift_factor)
        else:
            self.update_positions(self.environment.x_wind*wind_drift_factor,
                                  self.environment.y_wind*wind_drift_factor)

        ############################
        # Uncertainty / diffusion
        ############################
        if self.config['processes']['diffusion'] is True:
            # Current
            std_current_comp = self.config['drift']['current_uncertainty']
            if std_current_comp > 0:
                sigma_u = np.random.normal(0, std_current_comp,
                                           self.num_elements_active())
                sigma_v = np.random.normal(0, std_current_comp,
                                           self.num_elements_active())
            else:
                sigma_u = 0*self.environment.x_wind
                sigma_v = 0*self.environment.x_wind
            self.update_positions(sigma_u, sigma_v)

            # Wind
            std_wind_comp = self.config['drift']['wind_uncertainty']
            if wind_drift_factor > 0 and std_wind_comp > 0:
                sigma_u = np.random.normal(0, std_wind_comp*wind_drift_factor,
                                           self.num_elements_active())
                sigma_v = np.random.normal(0, std_wind_comp*wind_drift_factor,
                                           self.num_elements_active())
            else:
                sigma_u = 0*self.environment.x_wind
                sigma_v = 0*self.environment.x_wind
            self.update_positions(sigma_u, sigma_v)


    def vertical_advection(self, w_vel):
        """Move particles vertically according to vertical ocean current

            Vertical advection by ocean currents is small compared to termical velocity
            to do buoyancy and turbulent mixing
        """

    def vertical_mixing(self):
        """Mix particles vertically according to eddy diffusivity and buoyancy

            buoyancy is expressed as terminal velocity, which is the steady-state
            vertical velocity due to positive or negative buoyant behaviour.
            It is usually a function of particle density, diameter, and shape.
            
            Vertical particle displacemend du to turbulent mixing is calculated using 
            the "binned random walk scheme" (Thygessen and Aadlandsvik, 2007).
            The formulation of this scheme is copied from LADIM (IMR).
        """
        # if terminal_velocity is None:
        #     w = self.config['drift']['terminal_velocity']
        #else:
        #    w = terminal_velocity

        dz = 1. #m
        dt_mix = 2. #s

        # minimum height/maximum depth for each particle
        Zmin = -1.*self.environment.sea_floor_depth_below_sea_level

        # place particle in center of bin
        self.elements.z = np.round(self.elements.z/dz)*dz

        #avoid that elements are above surface / below bottom
        surface = np.where(self.elements.z > dz/2.)
        self.elements.z[surface] = dz/2.
        bottom = np.where(self.elements.z < Zmin)
        self.elements.z[bottom] = np.round(Zmin/dz)*dz + dz/2.

        # internal loop for fast time step of vertical mixing model
        # binned random walk needs faster time step compared to horizontal advection
        ntimes_mix = int(self.time_step.total_seconds()/dt_mix)
        logging.debug('Vertical mixing module:')
        logging.debug('turbulent diffusion with binned random walk scheme')
        logging.debug('using '+str(ntimes_mix)+' fast time steps of dt='+str(dt_mix)+'s')
        for i in range(0, ntimes_mix):
            # update terminal velocity according to environmental variables 
            self.update_terminal_velocity()
            w = self.elements.terminal_velocity
            K = self.environment.ocean_vertical_diffusivity

            #K1 = K # K at current depth
            #K2 = K # K at depth z-dz

            p = dt_mix * (2.0*K + dz*w)/(2.0*dz*dz)  #probability to rise
            q = dt_mix * (2.0*K - dz*w)/(2.0*dz*dz)  #probability to sink

            # check if probabilities are reasonable or wrong
            wrong = p+q > 1.
            if wrong.sum() > 0:
                logging.debug('WARNING! some elements have p+q>1.')

            RandKick = np.random.random(self.num_elements_active())
            up   = np.where(RandKick < p)
            down = np.where(RandKick > 1.0 - q)

            self.elements.z[up]   = self.elements.z[up]   + dz
            self.elements.z[down] = self.elements.z[down] - dz

            #avoid that elements are above surface / below bottom
            surface = np.where(self.elements.z > dz/2.)
            self.elements.z[surface] = dz/2.
            bottom = np.where(self.elements.z < Zmin)
            self.elements.z[bottom] = np.round(Zmin/dz)*dz + dz/2.

    def plot_vertical_distribution(self):
        """Function to plot vertical distribution of particles"""
        import matplotlib.pyplot as plt
        from matplotlib import dates
        dz=1.
        maxrange=-100
        hfmt = dates.DateFormatter('%d %b %Y %H:%M')
        fig = plt.figure()
        ax = fig.gca()
        #ax.xaxis.set_major_formatter(hfmt)
        #plt.xticks(rotation='vertical')
        #times = [self.start_time + n*self.time_step
        #         for n in range(self.steps + 1)]
        #data = self.history[prop].T[0:len(times), :]
        plt.hist(self.elements.z, bins=-maxrange/dz, range=[maxrange,0], orientation='horizontal')
        plt.title('Vertical particle distribution')
        plt.xlabel('number of particles')
        plt.ylabel('depth [m]')
        plt.subplots_adjust(bottom=.3)
        plt.grid()
        plt.show()


