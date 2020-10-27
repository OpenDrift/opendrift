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
from opendrift.models.oceandrift import OceanDrift
from opendrift.elements.passivetracer import PassiveTracer

# We add the property 'wind_drift_factor' to the element class
PassiveTracer.variables = PassiveTracer.add_variables([
    ('wind_drift_factor', {'dtype': np.float32,
                           'units': '1',
                           'default': 0.02}),
    ('terminal_velocity', {'dtype': np.float32,
                           'units': 'm/s',
                           'default': 0.01})])


class PlastDrift(OceanDrift):
    """Trajectory model based on the OpenDrift framework.

    Propagation of plastics particles with ocean currents and
    additional Stokes drift and wind drag.

    Developed at MET Norway.

    """

    ElementType = PassiveTracer

    required_variables = [
        'x_sea_water_velocity', 'y_sea_water_velocity',
        'sea_surface_wave_stokes_drift_x_velocity',
        'sea_surface_wave_stokes_drift_y_velocity',
        'sea_surface_wave_significant_height',
        'x_wind', 'y_wind',
        'ocean_vertical_diffusivity',
        'sea_floor_depth_below_sea_level']
    required_variables.append('land_binary_mask')

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_stokes_drift_x_velocity': 0,
                       'sea_surface_wave_stokes_drift_y_velocity': 0,
                       'x_wind': 0,
                       'y_wind': 0,
                       'sea_surface_wave_significant_height': 0,
                       'ocean_vertical_diffusivity': .02,
                       'sea_floor_depth_below_sea_level': 10000}

    required_profiles = ['ocean_vertical_diffusivity']


    def __init__(self, *args, **kwargs):

        # Call parent constructor
        super(PlastDrift, self).__init__(*args, **kwargs)

        # Add specific config settings
        self._add_config({
            # TODO: this option should be moved to OceanDrift
            'vertical_mixing:mixingmodel': {'type': 'enum',
                'enum': ['randomwalk', 'analytical'], 'default': 'analytical',
                'level': self.CONFIG_LEVEL_ADVANCED, 'description':
                    'Scheme to be used for vertical turbulent mixing'},
            })
    
        self._set_config_default('drift:vertical_advection', True)
        self._set_config_default('drift:use_tabularised_stokes_drift', True)
        self._set_config_default('general:coastline_action', 'previous')
        self._set_config_default('vertical_mixing:diffusivitymodel', 'windspeed_Large1994')

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


        if self.get_config('vertical_mixing:mixingmodel') == 'randomwalk':
            self.logger.debug('Turbulent mixing of particles using random walk')
            self.vertical_mixing()


        if self.get_config('vertical_mixing:mixingmodel') == 'analytical':
            self.logger.debug('Submerging according to wind')
            self.elements.z = -np.random.exponential(
                scale=self.environment.ocean_vertical_diffusivity/
                        self.elements.terminal_velocity,
                size=self.num_elements_active())
