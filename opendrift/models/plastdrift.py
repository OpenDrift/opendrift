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

import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED


class PlastElement(Lagrangian3DArray):
    variables = Lagrangian3DArray.add_variables([
        ('terminal_velocity', {'dtype': np.float32,
                               'units': 'm/s',
                               'level': CONFIG_LEVEL_ESSENTIAL,
            'description': 'Positive value means rising particles (positive buoyancy)',
                               'default': 0.01})])


class PlastDrift(OceanDrift):
    """Trajectory model based on the OpenDrift framework.

    Propagation of plastics particles with ocean currents and
    additional Stokes drift and wind drag.

    Developed at MET Norway.

    """

    ElementType = PlastElement

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_height': {'fallback': 0},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'ocean_vertical_diffusivity': {'fallback': 0.02, 'profiles': True},
        'ocean_mixed_layer_thickness': {'fallback': 50},
        'sea_floor_depth_below_sea_level': {'fallback': 10000},
        'land_binary_mask': {'fallback': None},
        }


    def __init__(self, *args, **kwargs):

        # Call parent constructor
        super(PlastDrift, self).__init__(*args, **kwargs)

        # Add specific config settings
        self._add_config({
            # TODO: this option should be moved to OceanDrift
            'vertical_mixing:mixingmodel': {'type': 'enum',
                'enum': ['randomwalk', 'analytical'], 'default': 'analytical',
                'level': CONFIG_LEVEL_ADVANCED, 'description':
                    'Scheme to be used for vertical turbulent mixing'},
            })

        self._set_config_default('drift:vertical_mixing', True)
        self._set_config_default('drift:vertical_advection', True)
        self._set_config_default('drift:use_tabularised_stokes_drift', True)
        self._set_config_default('general:coastline_action', 'previous')
        self._set_config_default('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983')

    def update(self):
        """Update positions and properties of elements."""

        # Simply move particles with ambient current
        self.advect_ocean_current()

        self.update_particle_depth()

        # Advect particles due to Stokes drift
        self.stokes_drift()

        # Advect particles due to wind-induced shear near surface
        self.advect_wind()

    def update_particle_depth(self):

        if self.get_config('drift:vertical_mixing') is True:

            if self.get_config('vertical_mixing:mixingmodel') == 'randomwalk':
                logger.debug('Turbulent mixing of particles using random walk')
                self.vertical_mixing()

            if self.get_config('vertical_mixing:mixingmodel') == 'analytical':
                logger.debug('Submerging according to wind')
                self.elements.z = -np.random.exponential(
                    scale=self.environment.ocean_vertical_diffusivity/
                            self.elements.terminal_velocity,
                    size=self.num_elements_active())
