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

from openoil import OpenOil, Oil
from opendrift3D import OpenDrift3DSimulation, Lagrangian3DArray


# Defining the oil element properties
class Oil3D(Oil):
    """Extending Oil class with variables relevant for the vertical."""

    variables = Oil.add_variables([
        ('diameter', {'dtype': np.float32,
                       'units': 'm',
                       'default': 0.0001})
        ])


class OpenOil3D(OpenDrift3DSimulation, OpenOil):  # Multiple inheritance
    """Open source oil trajectory model based on the OpenDrift framework.

        Developed at MET Norway based on oil weathering parameterisations
        found in open/published litterature.

        Under construction.
    """

    ElementType = Oil3D

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


    # Read oil types from file (presently only for illustrative effect)
    oil_types = str([str(l.strip()) for l in open(
                    os.path.dirname(os.path.realpath(__file__))
                    + '/oil_types.txt').readlines()])[1:-1]
    default_oil = oil_types.split(',')[0].strip()


    # Configuration

    configspec = '''
        [input]
            readers = list(min=1, default=list(''))
            [[spill]]
                longitude = float(min=-360, max=360, default=5)
                latitude = float(min=-90, max=90, default=60)
                time = string(default='%s')
                oil_type = option(%s, default=%s)
        [processes]
            turbulentmixing = boolean(default=True)
            verticaladvection = boolean(default=False)
            dispersion = boolean(default=True)
            diffusion = boolean(default=True)
            evaporation = boolean(default=True)
            emulsification = boolean(default=True)
        [drift]
            wind_drift_factor = float(min=0, max=1, default=0.02)
            current_uncertainty = float(min=0, max=5, default=.1)
            wind_uncertainty = float(min=0, max=5, default=1)
            relative_wind = boolean(default=True)
        [turbulentmixing]
            timestep = float(min=0.1, max=3600, default=1.)
            verticalresolution = float(min=0.01, max=10, default = 1.)
            diffusivitymodel = string(default='environment')

    ''' % (datetime.now().strftime('%Y-%d-%m %H:00'), oil_types, default_oil)

    def __init__(self, *args, **kwargs):

        # Read oil properties from file
        self.oiltype_file = os.path.dirname(os.path.realpath(__file__)) + \
                                '/oilprop.dat'
        oilprop = open(self.oiltype_file)
        oiltypes = []
        linenumbers = []
        for i, line in enumerate(oilprop.readlines()):
            if line[0].isalpha():
                oiltype = line.strip()[:-2].strip()
                oiltypes.append(oiltype)
                linenumbers.append(i)
        oiltypes, linenumbers = zip(*sorted(zip(oiltypes, linenumbers)))
        self.oiltypes = oiltypes
        self.oiltypes_linenumbers = linenumbers

        # Calling general constructor of parent class
        super(OpenOil3D, self).__init__(*args, **kwargs)

    def update_terminal_velocity(self):
        """Calculate terminal velocity for oil droplets

        according to 
        Tkalich et al. (2002): Vertical mixing of oil droplets
                               by breaking waves
        Marine Pollution Bulletin 44, 1219-1229
        """
        g = 9.81 # ms-2

        r = self.elements.diameter # particle radius
        T0 = self.environment.sea_water_temperature
        S0 = self.environment.sea_water_salinity
        rho_oil = self.elements.density
        rho_water = self.sea_water_density(T=T0, S=S0)

        # dynamic water viscosity
        my_w = 0.001*(1.7915 - 0.0538*T0 + 0.007*(T0**(2.0)) - 0.0023*S0) # ~0.0014 kg m-1 s-1
        # kinemativ water viscosity
        ny_w = my_w / rho_water
        # 
        rhopr = rho_oil/rho_water

        # terminal velocity for low Reynolds numbers
        kw = 2*g*(1-rhopr)/(9*ny_w)
        W = kw * r**2
        #print r[0:10], 'r'
        #print ny_w[0:10], 'ny_w'
        #print W[0:10], 'W before'

        #check if we are in a high Reynolds number regime
        Re = 2*r*W/ny_w
        highRe = np.where(Re > 50) 
        
        # Terminal velocity in high Reynolds numbers
        kw = (16*g*(1-rhopr)/3)**0.5
        W2 = kw*r**0.5

        W[highRe] = W2[highRe]
        #print W[0:10], 'W after'
        self.elements.terminal_velocity = W 

    def update(self):
        """Update positions and properties of oil particles."""

        self.elements.age_seconds += self.time_step.total_seconds()

        # Calculate windspeed, which is needed for emulsification,
        # and dispersion (if wave height is not given)
        windspeed = np.sqrt(self.environment.x_wind**2 +
                            self.environment.y_wind**2)

        ################
        ## Evaporation
        ################
        reference_wind = 10  # Hardcoded for now, should be read from file
        if self.config['processes']['evaporation'] is True:
            # Store evaporation fraction at beginning of timestep
            fraction_evaporated_previous = self.elements.fraction_evaporated

            # Evaporate only elements at surface
            at_surface = (self.elements.z == 0)
            if np.isscalar(at_surface):
                at_surface = at_surface*np.ones(self.num_elements_active(),
                                                dtype=bool)
            Urel = windspeed/reference_wind  # Relative wind
            h = 2  # Film thickness in mm, harcoded for now

            # Calculate exposure time
            #   presently without compensation for sea temperature
            delta_exposure_seconds = \
                (self.model.reference_thickness/h)*Urel * \
                self.time_step.total_seconds()
            if np.isscalar(self.elements.age_exposure_seconds):
                self.elements.age_exposure_seconds += delta_exposure_seconds
            else:
                self.elements.age_exposure_seconds[at_surface] += \
                    delta_exposure_seconds[at_surface]

            self.elements.fraction_evaporated = np.interp(
                self.elements.age_exposure_seconds,
                self.model.tref, self.model.fref)
            self.mass_evaporated = \
                self.elements.mass_oil*self.elements.fraction_evaporated
            # Remove evaporated part from mass_oil
            mass_evaporated_timestep = self.elements.mass_oil*(
                self.elements.fraction_evaporated -
                fraction_evaporated_previous)
            self.elements.mass_oil -= mass_evaporated_timestep
            self.elements.mass_evaporated += mass_evaporated_timestep

            # Evaporation probability equals the difference in fraction
            # evaporated at this timestep compared to previous timestep,
            # divided by the remaining fraction of the particle at
            # previous timestep
            evaporation_probability = ((self.elements.fraction_evaporated -
                                        fraction_evaporated_previous) /
                                       (1 - fraction_evaporated_previous))
            evaporation_probability[~at_surface] = 0
            evaporated_indices = (evaporation_probability >
                                  np.random.rand(self.num_elements_active(),))
            self.deactivate_elements(evaporated_indices, reason='evaporated')

        ##################
        # Emulsification
        ##################
        if self.config['processes']['emulsification'] is True:
            # Apparent emulsion age of particles
            Urel = windspeed/reference_wind  # Relative wind
            self.elements.age_emulsion_seconds += \
                Urel*self.time_step.total_seconds()

            self.elements.water_content = np.interp(
                self.elements.age_emulsion_seconds,
                self.model.tref, self.model.wmax)

        ###############
        # Dispersion
        ###############
        if self.config['processes']['dispersion'] is True:

            # From NOAA PyGnome model:
            # https://github.com/NOAA-ORR-ERD/PyGnome/
            v_entrain = 3.9E-8
            sea_water_density = 1028
            fraction_breaking_waves = 0.02
            wave_significant_height = \
                self.environment.sea_surface_wave_significant_height
            wave_significant_height[wave_significant_height == 0] = \
                0.0246*windspeed[wave_significant_height == 0]**2
            dissipation_wave_energy = \
                (0.0034*sea_water_density*9.81*wave_significant_height**2)
            c_disp = np.power(dissipation_wave_energy, 0.57) * \
                fraction_breaking_waves
            # Roy's constant
            C_Roy = 2400.0 * np.exp(-73.682*np.sqrt(
                self.elements.viscosity/self.elements.density))

            q_disp = C_Roy * c_disp * v_entrain / self.elements.density

            oil_mass_loss = (q_disp * self.time_step.total_seconds() *
                             self.elements.density)*self.elements.mass_oil

            self.elements.mass_oil -= oil_mass_loss
            self.elements.mass_dispersed += oil_mass_loss

            #self.elements.z = \
            #    self.environment.sea_surface_wave_significant_height

            ## Marks R. (1987), Marine aerosols and whitecaps in the
            ## North Atlantic and Greenland sea regions
            ## Deutsche Hydrografische Zeitschrift, Vol 40, Issue 2 , pp 71-79
            #whitecap_coverage = (2.54E-4)*np.power(windspeed, 3.58)

            ## Martinsen et al. (1994), The operational
            ## oil drift system at DNMI
            ## DNMI Technical report No 125, 51
            #wave_period = 3.85*np.sqrt(
            #    self.environment.sea_surface_wave_significant_height)  # sec

            #time_between_breaking_events = 3.85/whitecap_coverage  # TBC

            #rho_w = 1025  # kg/m3
            #wave_period[wave_period == 0] = 5  # NB: temporal fix for no waves
            #dissipation_energy = 0.0034*rho_w*9.81*(wave_period**2)
            #dsize_coeff = 2100  # Wettre et al., Appendix A
            #c_oil = dsize_coeff*np.power(self.elements.viscosity, -0.4)
            ## Random numbers between 0 and 1:
            #p = np.random.rand(self.num_elements_active(), )
            #oil_per_unit_surface = 1
            #droplet_size = np.power((p*oil_per_unit_surface)/(c_oil *
            #                        np.power(dissipation_energy, 0.57)),
            #                        1/1.17)
            #self.deactivate_elements(droplet_size < 5E-7, reason='dispersed')

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
        # Deactivate elements hitting coast
        ####################################
        self.deactivate_elements(self.environment.land_binary_mask == 1,
                                 reason='stranded')

        ##########################################
        # Do not let particles go below seafloor
        ##########################################
        self.elements.z[np.where(self.elements.z < \
            -self.environment.sea_floor_depth_below_sea_level)] = \
            -self.environment.sea_floor_depth_below_sea_level

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
