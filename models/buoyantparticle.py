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
from datetime import datetime

from opendrift import OpenDriftSimulation
from elements import LagrangianArray


# Defining the oil element properties
class BuoyantParticle(LagrangianArray):
    """Extending LagrangianArray to a passive particle subject to mixing that can be positively or negatively buoyant"""

    variables = LagrangianArray.add_variables([
        ('buoyancy', {'dtype': np.float32,
                      'units': '',
                      'default': 0.}),
        ('diameter', {'dtype': np.float32,
                       'units': 'm',
                       'default': 0.001}),
        ('density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': 800}),
        ('terminal_velocity', {'dtype': np.float32,
                     'units': 'm/s',
                     'default': 0.}),
        ('age_seconds', {'dtype': np.float32,
                         'units': 's',
                         'default': 0}),])
    def update_terminal_velocity(self):
        """Calculate terminal velocity due to bouyancy from own properties
        and environmental variables"""
        self.terminal_velocity = 0.001


class OpenBuoyantParticle(OpenDriftSimulation):
    """Open source buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway
        
        Generic module for particles that are subject to vertical turbulent mixing
        with the possibility for positive of negative buoyancy

        Particles could be pollution (e.g. oil droplets), plankton, or sediments

        Under construction.
    """

    ElementType = BuoyantParticle

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'sea_surface_wave_significant_height',
                          'sea_surface_wave_to_direction',
                          'sea_ice_area_fraction',
                          'x_wind', 'y_wind', 'land_binary_mask',
                          'ocean_vertical_diffusivity'
                          #'upward_sea_water_velocity'
                         ]

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       'sea_surface_wave_to_direction': 0,
                       'sea_ice_area_fraction': 0,
                       'x_wind': 0, 'y_wind': 0,
                       'ocean_vertical_diffusivity': 0.05 #m2s-1
                       #'upward_sea_water_velocity': 0
                      }

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'missing_data': 'gray', 'stranded': 'red',
                     'evaporated': 'yellow', 'dispersed': 'magenta'}

    # Read oil types from file (presently only for illustrative effect)
    #oil_types = str([str(l.strip()) for l in open(
    #                os.path.dirname(os.path.realpath(__file__))
    #                + '/oil_types.txt').readlines()])[1:-1]
    #default_oil = oil_types.split(',')[0].strip()


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

        # Read oil properties from file
        #self.oiltype_file = os.path.dirname(os.path.realpath(__file__)) + \
        #                        '/oilprop.dat'
        #oilprop = open(self.oiltype_file)
        #oiltypes = []
        #linenumbers = []
        #for i, line in enumerate(oilprop.readlines()):
        #    if line[0].isalpha():
        #        oiltype = line.strip()[:-2].strip()
        #        oiltypes.append(oiltype)
        #        linenumbers.append(i)
        #oiltypes, linenumbers = zip(*sorted(zip(oiltypes, linenumbers)))
        #self.oiltypes = oiltypes
        #self.oiltypes_linenumbers = linenumbers

        # Calling general constructor of parent class
        super(OpenBuoyantParticle, self).__init__(*args, **kwargs)

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


            #sea_water_density = 1028

            # terminal velocity is the uprising/sedimentation velocity due to buoyancy
            # it can be calculated from density and diameter
            self.elements.update_terminal_velocity() # 0.001m/s (1mm/s)

            #self.vertical_mixing(self.environment.ocean_vertical_diffusivity, self.elements.terminal_velocity)
            self.vertical_mixing(self.environment.ocean_vertical_diffusivity, 0.001)



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

        ########################################
        ## Deactivate elements hitting sea ice
        ########################################
        self.deactivate_elements(self.environment.sea_ice_area_fraction > 0.6,
                                 reason='oil-in-ice')

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

    def plot_vertical_distribution(self):
        #def plot_oil_budget(self):
        import matplotlib.pyplot as plt
        mass_oil, status = self.get_property('mass_oil')
        if 'stranded' not in self.status_categories:
            self.status_categories.append('stranded')
        mass_active = np.ma.sum(np.ma.masked_where(
                status==self.status_categories.index('stranded'),
                    mass_oil), axis=1)
        mass_stranded = np.ma.sum(np.ma.masked_where(
            status!=self.status_categories.index('stranded'),
            mass_oil), axis=1)
        mass_evaporated, status = self.get_property('mass_evaporated')
        mass_evaporated = np.sum(mass_evaporated, axis=1)
        mass_dispersed, status = self.get_property('mass_dispersed')
        mass_dispersed = np.sum(mass_dispersed, axis=1)

        budget = np.row_stack((mass_dispersed, mass_active,
                               mass_stranded, mass_evaporated))
        budget = np.cumsum(budget, axis=0)

        time, time_relative = self.get_time_array()
        time = np.array([t.total_seconds()/3600. for t in time_relative])
        fig = plt.figure()
        # Left axis showing oil mass
        ax1 = fig.add_subplot(111)
        # Hack: make some emply plots since fill_between does not support label
        ax1.add_patch(plt.Rectangle((0, 0), 0, 0,
                      color='salmon', label='dispersed'))
        ax1.add_patch(plt.Rectangle((0, 0), 0, 0,
                      color='royalblue', label='in water'))
        ax1.add_patch(plt.Rectangle((0, 0), 0, 0,
                      color='black', label='stranded'))
        ax1.add_patch(plt.Rectangle((0, 0), 0, 0,
                      color='skyblue', label='evaporated'))

        ax1.fill_between(time, 0, budget[0,:], facecolor='salmon')
        ax1.fill_between(time, budget[0,:], budget[1,:],
                         facecolor='royalblue')
        ax1.fill_between(time, budget[1,:], budget[2,:],
                         facecolor='black')
        ax1.fill_between(time, budget[2,:], budget[3,:],
                         facecolor='skyblue')
        ax1.set_ylim([0,budget.max()])
        ax1.set_xlim([0,time.max()])
        ax1.set_ylabel('Mass oil  [%s]' %
                        self.elements.variables['mass_oil']['units'])
        ax1.set_xlabel('Time  [hours]')
        # Right axis showing percent
        ax2 = ax1.twinx()
        ax2.set_ylim([0,100])
        ax2.set_ylabel('Percent')
        plt.title('%s - %s to %s' %
                  (self.model.oiltype,
                   self.start_time.strftime('%Y-%m-%d %H:%M'),
                   self.time.strftime('%Y-%m-%d %H:%M')))
        # Shrink current axis's height by 10% on the bottom
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        ax2.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        ax1.legend(bbox_to_anchor=(0., -0.10, 1., -0.03), loc=1,
                   ncol=4, mode="expand", borderaxespad=0.)
        plt.show()

        self.schedule_elements(self.ElementType(**kwargs), oil_time)
