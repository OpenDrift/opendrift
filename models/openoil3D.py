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

    required_profiles = ['sea_water_temperature',
                         'ocean_vertical_diffusivity']
    required_profiles_z_range = [-50, 0]  # The depth range (in m) which
                                          # profiles shall cover

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       'sea_surface_wave_to_direction': 0,
                       'sea_ice_area_fraction': 0,
                       'x_wind': 0, 'y_wind': 0,
                       'sea_floor_depth_below_sea_level': 100,
                       'ocean_vertical_diffusivity': 0.02,  # m2s-1
                       'sea_water_temperature': 10.,
                       'sea_water_salinity': 34.
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
            scheme = option('euler', 'runge-kutta', default='euler')
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
        g = 9.81  # ms-2

        r = self.elements.diameter  # particle radius
        T0 = self.environment.sea_water_temperature
        S0 = self.environment.sea_water_salinity
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

        # Oil weathering (inherited from OpenOil)
        self.oil_weathering()

        # Turbulent Mixing
        if self.config['processes']['turbulentmixing'] is True:
            self.update_terminal_velocity()
            self.vertical_mixing()

        # Vertical advection
        if self.config['processes']['verticaladvection'] is True:

            self.vertical_advection(self.environment.upward_sea_water_velocity)

        # Horizontal advection (inherited from OpenOil)
        self.advect_oil()
