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

from opendrift3D import OpenDrift3DSimulation, Lagrangian3DArray
from elements import LagrangianArray


# Defining the oil element properties
class PelagicEgg(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for pelagic eggs
    """

    variables = LagrangianArray.add_variables([
        ('diameter', {'dtype': np.float32,
                       'units': 'm',
                       'default': 0.0014}), # for NEA Cod
        ('neutral_buoyancy_salinity', {'dtype': np.float32,
                       'units': '[]',
                       'default': 31.25}), # for NEA Cod
        ('density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': 1028.}),
        ('age_seconds', {'dtype': np.float32,
                         'units': 's',
                         'default': 0.}),        
        ('hatched', {'dtype': np.float32,
                     'units': '',
                     'default': 0.})])

class PelagicEggDrift(OpenDrift3DSimulation):
    """Open source buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway
        
        Generic module for particles that are subject to vertical turbulent mixing
        with the possibility for positive of negative buoyancy

        Particles could be pollution (e.g. oil droplets), plankton, or sediments

        Under construction.
    """

    ElementType = PelagicEgg

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'sea_surface_wave_significant_height',
                          'sea_surface_wave_to_direction',
                          'sea_ice_area_fraction',
                          'x_wind', 'y_wind', 'land_binary_mask',
                          'sea_floor_depth_below_sea_level',
                          'ocean_vertical_diffusivity',
                          'sea_water_temperature',
                          'sea_water_salinity'
                          #'upward_sea_water_velocity'
                         ]

    # Vertical profiles of the following parameters will be available in
    # dictionary self.environment.vertical_profiles
    # E.g. self.environment_profiles['x_sea_water_velocity']
    # will be an array of size [vertical_levels, num_elements]
    # The vertical levels are available as
    # self.environment_profiles['z'] or
    # self.environment_profiles['sigma'] (not yet implemented) 
    required_profiles = ['sea_water_temperature', 'sea_water_salinity']
    required_profiles_z_range = [-50, 0]  # The depth range (in m) which
                                          # profiles shall cover

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       'sea_surface_wave_to_direction': 0,
                       'sea_ice_area_fraction': 0,
                       'x_wind': 0, 'y_wind': 0,
                       'sea_floor_depth_below_sea_level': 100,
                       'ocean_vertical_diffusivity': 0.02, #m2s-1
                       'sea_water_temperature': 10.,
                       'sea_water_salinity' : 34.
                       #'upward_sea_water_velocity': 0
                      }

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'hatched': 'red', 'eaten': 'yellow', 'died': 'magenta'}

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
        [turbulentmixing]
            timestep = float(min=0.1, max=3600, default=1.)
            verticalresolution = float(min=0.01, max=10, default = 1.)
            diffusivitymodel = string(default='environment') 

    ''' % (datetime.now().strftime('%Y-%d-%m %H:00'))

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(PelagicEggDrift, self).__init__(*args, **kwargs)

#    def update_terminal_velocity(self):
#        """Calculate terminal velocity due to bouyancy from own properties
#        and environmental variables"""
#        self.elements.terminal_velocity = 0.005
        
    def update_terminal_velocity(self):
        """Calculate terminal velocity for Pelagic Egg

        according to 
        S. Sundby (1983): A one-dimensional model for the vertical distribution
        of pelagic fish eggs in the mixed layer
        Deep Sea Research (30) pp. 645-661

        Method copied from ibm.f90 module of LADIM:
        Vikebo, F., S. Sundby, B. Aadlandsvik and O. Otteraa (2007),
        Fish. Oceanogr. (16) pp. 216-228
        """
        g = 9.81 # ms-2

        # Pelagic Egg properties that determine buoyancy
        eggsize = self.elements.diameter # 0.0014 for NEA Cod
        eggsalinity = self.elements.neutral_buoyancy_salinity # 31.25 for NEA Cod

        T0 = self.environment.sea_water_temperature
        S0 = self.environment.sea_water_salinity

        # The density difference bettwen a pelagic egg and the ambient water
        # is regulated by their salinity difference through the equation of state for sea water.
        # The Egg has the same temperature as the abient water and its salinity is
        # regulated by osmosis through the egg shell.
        DENSw   = sea_water_density(T=T0, S=S0)
        DENSegg = sea_water_density(T=T0, S=eggsalinity)
        dr = DENSw-DENSegg # density difference

        # water viscosity
        my_w = 0.001*(1.7915 - 0.0538*T0 + 0.007*(T0**(2.0)) - 0.0023*S0) # ~0.0014 kg m-1 s-1

        # terminal velocity for low Reynolds numbers
        W = (1.0/my_w)*(1.0/18.0)*g*eggsize**2 * dr

        #check if we are in a Reynolds regime where Re > 0.5
        highRe = np.where(W*1000*eggsize/my_w > 0.5) 
        
        # Use empirical equations for terminal velocity in high Reynolds numbers
        # empirical equations have length units in cm!
        my_w = 0.01854*np.exp(-0.02783*T0) #in cm2/s 
        d0 = (eggsize*100) - 0.4*(9.0*my_w**2 /(100*g) *DENSw/dr)**(1.0/3.0) #cm 
        W2 = 19.0* d0* (0.001*dr)**(2.0/3.0)*(my_w*0.001*DENSw)**(-1.0/3.0) #cm/s
        W2 = W2/100. #back to m/s

        W[highRe] = W2[highRe]
        self.elements.terminal_velocity = W

    def update(self):
        """Update positions and properties of buoyant particles."""

        self.elements.age_seconds += self.time_step.total_seconds()

        # Calculate windspeed, which is needed for emulsification,
        # and dispersion (if wave height is not given)
        # wind speed can also be used to parameterize eddy diffusivity
        windspeed = np.sqrt(self.environment.x_wind**2 +
                            self.environment.y_wind**2)

        # For illustration: plot vertical profile for a single element
        if self.steps == 5 and self.environment_profiles is not None:
            # We plot profile at a single timestep only
            import matplotlib.pyplot as plt
            param = self.required_profiles[0]  # take the first one
            elem_num = 0
            plt.plot(self.environment_profiles[param][:,elem_num],
                     self.environment_profiles['z'], '-*')
            plt.plot(self.environment[param][elem_num],
                     self.elements.z[elem_num], 'r*', markersize=15)
            plt.ylabel('Depth  [m]')
            plt.xlabel(param)
            plt.title('%s  (%fN, %fE)' %
                        (self.time, self.elements.lat[elem_num],
                         self.elements.lon[elem_num]))
            plt.gca().set_yticks(self.environment_profiles['z'])
            plt.grid('on')
            plt.show()

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

def sea_water_density(T=10.,S=35.):
    '''The function gives the density of seawater at one atmosphere
    pressure as given in :

    N.P. Fofonoff and R.C. Millard Jr.,1983,
    Unesco technical papers in marine science no. 44.
    
    S   = Salinity in promille of the seawater
    T   = Temperature of the seawater in degrees Celsius
    '''

    R4=4.8314E-04
    DR350=28.106331

    #Pure water density at atmospheric pressure
    #Bigg P.H. (1967) BR. J. Applied Physics pp.:521-537

    R1 = ((((6.536332E-09*T-1.120083E-06)*T+1.001685E-04)*T-9.095290E-03)*T+6.793952E-02)*T-28.263737
    
    #Seawater density at atmospheric pressure
    #coefficients involving salinity :

    R2 = (((5.3875E-09*T-8.2467E-07)*T+7.6438E-05)*T-4.0899E-03)*T +8.24493E-01

    R3 = (-1.6546E-06*T+1.0227E-04)*T-5.72466E-03

    #International one-atmosphere equation of state of seawater :

    SIG = R1 + (R4*S + R3*np.sqrt(S) + R2)*S
    Dens0 = SIG + DR350 + 1000.
    return Dens0
