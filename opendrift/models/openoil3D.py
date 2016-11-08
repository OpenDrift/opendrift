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

from opendrift.models.openoil import OpenOil, Oil
from opendrift.models.opendrift3D import OpenDrift3DSimulation


# Defining the oil element properties
class Oil3D(Oil):
    """Extending Oil class with variables relevant for the vertical."""

    variables = Oil.add_variables([
        # Entrainment length scale, see Tkalich and Chan (2002)
        ('entrainment_length_scale', {'dtype': np.float32,
                                      'units': 'm',
                                      'default': 0.1}),
        ('diameter', {'dtype': np.float32,  # Particle diameter
                      'units': 'm',
                      'default': 0.00001})
        ])


class OpenOil3D(OpenDrift3DSimulation, OpenOil):  # Multiple inheritance
    """Open source oil trajectory model based on the OpenDrift framework.

        Developed at MET Norway based on oil weathering parameterisations
        found in open/published litterature.

        Under construction.
    """

    ElementType = Oil3D

    required_variables = [
        'x_sea_water_velocity', 'y_sea_water_velocity',
        'sea_surface_wave_significant_height',
        'sea_surface_wave_to_direction',
        'sea_surface_wave_stokes_drift_x_velocity',
        'sea_surface_wave_stokes_drift_y_velocity',
        'sea_surface_wave_period_at_variance_spectral_density_maximum',
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
        'sea_ice_area_fraction',
        'x_wind', 'y_wind', 'land_binary_mask',
        'sea_floor_depth_below_sea_level',
        'ocean_vertical_diffusivity',
        'sea_water_temperature',
        'sea_water_salinity',
        'upward_sea_water_velocity'
        ]

    required_profiles = ['sea_water_temperature',
                         'sea_water_salinity',
                         'ocean_vertical_diffusivity']
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-120, 0]

    fallback_values = {
        'x_sea_water_velocity': 0,
        'y_sea_water_velocity': 0,
        'sea_surface_wave_significant_height': 0,
        'sea_surface_wave_to_direction': 0,
        'sea_surface_wave_stokes_drift_x_velocity': 0,
        'sea_surface_wave_stokes_drift_y_velocity': 0,
        'sea_surface_wave_period_at_variance_spectral_density_maximum': 0,
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 0,
        'sea_ice_area_fraction': 0,
        'x_wind': 0, 'y_wind': 0,
        'sea_floor_depth_below_sea_level': 10000,
        'ocean_vertical_diffusivity': 0.02,  # m2s-1
        'sea_water_temperature': 10.,
        'sea_water_salinity': 34.,
        'upward_sea_water_velocity': 0
        }

    # Read oil types from file (presently only for illustrative effect)
    oil_types = str([str(l.strip()) for l in open(
                    os.path.dirname(os.path.realpath(__file__)) +
                    '/oil_types.txt').readlines()])[1:-1]
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
            scheme = option('euler', 'runge-kutta', default='euler')
            wind_drift_factor = float(min=0, max=1, default=0.02)
            current_uncertainty = float(min=0, max=5, default=.1)
            wind_uncertainty = float(min=0, max=5, default=1)
            relative_wind = boolean(default=True)
        [turbulentmixing]
            timestep = float(min=0.1, max=3600, default=1.)
            verticalresolution = float(min=0.01, max=10, default = 1.)
            diffusivitymodel = string(default='environment')
            TSprofiles = boolean(default=False)

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

    def particle_radius(self):
        """Calculate radius of entained particles.

        Per now a fixed radius, should later use a distribution.
        """

        # Delvigne and Sweeney (1988)
        # rmax = 1818*np.power(self.wave_energy_dissipation(), -0.5)* \
        #             np.power(self.elements.viscosity, 0.34) / 1000000

        # r = np.random.uniform(0, rmax, self.num_elements_active())
        return self.elements.diameter/2  # Hardcoded diameter

    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        """Calculate terminal velocity for oil droplets

        according to
        Tkalich et al. (2002): Vertical mixing of oil droplets
                               by breaking waves
        Marine Pollution Bulletin 44, 1219-1229

        If profiles of temperature and salt are passed into this function,
        they will be interpolated from the profiles.
        if not, T,S will be fetched from reader.
        """
        g = 9.81  # ms-2

        r = self.particle_radius()*2

        # prepare interpolation of temp, salt

        if not (Tprofiles is None and Sprofiles is None):
            if z_index is None:
                z_i = range(Tprofiles.shape[0])  # evtl. move out of loop
                # evtl. move out of loop
                z_index = interp1d(-self.environment_profiles['z'],
                                   z_i, bounds_error=False)
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.int), 0)
            lower = np.minimum(upper+1, Tprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)

        # do interpolation of temp, salt if profiles were passed into
        # this function, if not, use reader by calling self.environment
        if Tprofiles is None:
            T0 = self.environment.sea_water_temperature
        else:
            T0 = Tprofiles[upper, range(Tprofiles.shape[1])] * \
                weight_upper + \
                Tprofiles[lower, range(Tprofiles.shape[1])] * \
                (1-weight_upper)
        if Sprofiles is None:
            S0 = self.environment.sea_water_salinity
        else:
            S0 = Sprofiles[upper, range(Sprofiles.shape[1])] * \
                weight_upper + \
                Sprofiles[lower, range(Sprofiles.shape[1])] * \
                (1-weight_upper)

        rho_oil = self.elements.density
        rho_water = self.sea_water_density(T=T0, S=S0)

        # dynamic water viscosity
        my_w = 0.001*(1.7915 - 0.0538*T0 + 0.007*(T0**(2.0)) - 0.0023*S0)
        # ~0.0014 kg m-1 s-1
        # kinemativ water viscosity
        ny_w = my_w / rho_water
        rhopr = rho_oil/rho_water

        # terminal velocity for low Reynolds numbers
        kw = 2*g*(1-rhopr)/(9*ny_w)
        W = kw * r**2

        # check if we are in a high Reynolds number regime
        Re = 2*r*W/ny_w
        highRe = np.where(Re > 50)

        # Terminal velocity in high Reynolds numbers
        kw = (16*g*(1-rhopr)/3)**0.5
        W2 = kw*r**0.5

        W[highRe] = W2[highRe]
        self.elements.terminal_velocity = W

    def oil_wave_entrainment_rate(self):
        kb = 0.4
        omega = (2.*np.pi)/self.wave_period()
        gamma = self.wave_damping_coefficient()
        alpha = 1.5
        Low = self.elements.entrainment_length_scale
        entrainment_rate = \
            kb*omega*gamma*self.significant_wave_height() / \
            (16*alpha*Low)
        return entrainment_rate

    def prepare_vertical_mixing(self):
        '''Calculate entrainment probability before main loop'''
        self.oil_entrainment_probability = \
            self.oil_wave_entrainment_rate()*\
                self.config['turbulentmixing']['timestep']

    def surface_interaction(self, time_step_seconds, alpha=1.5):
        """Mix surface oil into water column."""

        # Place particles above surface to exactly 0
        surface = self.elements.z >= 0
        self.elements.z[surface] = 0

        # Entrain oil into uppermost layer (whitecapping from waves)
        # TODO: optimise this by only calculate for surface elements
        random_number = np.random.uniform(0, 1, len(self.elements.z))
        entrained = np.logical_and(surface,
                        random_number<self.oil_entrainment_probability)
        self.elements.z[entrained] = \
            -self.config['turbulentmixing']['verticalresolution']/2.

    def resurface_elements(self, minimum_depth=None):
        """Oil elements reaching surface (or above) form slick, not droplet"""
        surface = np.where(self.elements.z >= 0)[0]
        self.elements.z[surface] = 0

    def update(self):
        """Update positions and properties of oil particles."""

        # Oil weathering (inherited from OpenOil)
        self.oil_weathering()

        # Turbulent Mixing
        if self.config['processes']['turbulentmixing'] is True:
            self.update_terminal_velocity()
            self.vertical_mixing()

        # Vertical advection
        if self.config['processes']['verticaladvection'] is True:
            self.vertical_advection()

        # Horizontal advection (inherited from OpenOil)
        self.advect_oil()
