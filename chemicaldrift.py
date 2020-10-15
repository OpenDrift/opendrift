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
# Copyright 2020, Knut-Frode Dagestad, MET Norway

"""
ChemicalDrift is an OpenDrift module for drift and settling of chemicals.
Based on work by Simon Weppe, MetOcean Solutions Ltd.
"""

import numpy as np
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.oceandrift import Lagrangian3DArray

class ChemicalElement(Lagrangian3DArray):
    variables = Lagrangian3DArray.add_variables([
        ('settled', {'dtype': np.int16,  # 0 is active, 1 is settled
                     'units': '1',
                     'default': 0}),
        ('terminal_velocity', {'dtype': np.float32,
                               'units': 'm/s',
                               'default': -0.001})  # 1 mm/s negative buoyancy
        ])


class ChemicalDrift(OceanDrift):
    """Model for chemical drift, under development
    """

    ElementType = ChemicalElement

    required_variables = [
            'x_sea_water_velocity', 'y_sea_water_velocity',
            'upward_sea_water_velocity',
            'x_wind', 'y_wind',
            'sea_surface_wave_stokes_drift_x_velocity',
            'sea_surface_wave_stokes_drift_y_velocity',
            'sea_surface_wave_period_at_variance_spectral_density_maximum',
            'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
            'land_binary_mask',
            'ocean_vertical_diffusivity',
            'sea_floor_depth_below_sea_level'
            ]

    fallback_values = {
            'x_sea_water_velocity': 0, 'y_sea_water_velocity': 0,
            'upward_sea_water_velocity': 0,
            'x_wind': 0, 'y_wind': 0,
            'sea_surface_wave_stokes_drift_x_velocity': 0,
            'sea_surface_wave_stokes_drift_y_velocity': 0,
            'sea_surface_wave_period_at_variance_spectral_density_maximum': 0,
            'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 0,
            'ocean_vertical_diffusivity': .02,
            'sea_floor_depth_below_sea_level': 10000
            }

    def __init__(self, *args, **kwargs):
        """ Constructor of ChemicalDrift module
        """

        super(ChemicalDrift, self).__init__(*args, **kwargs)

        # By default, chemicals do not strand towards coastline
        # TODO: A more sophisticated stranding algorithm is needed
        self.set_config('general:coastline_action', 'previous')

        # Vertical mixing is enabled as default
        self.set_config('drift:vertical_mixing', True)

    def update(self):
        """Update positions and properties of chemical particles.
        """

        # Advecting here all elements, but want to soon add 
        # possibility of not moving settled elements, until
        # they are resuspended. May then need to send a boolean
        # array to advection methods below
        self.advect_ocean_current()

        self.vertical_advection()

        self.advect_wind()  # Wind shear in upper 10cm of ocean

        self.stokes_drift()

        self.vertical_mixing()  # Including buoyancy and settling

        self.resuspension()

    def bottom_interaction(self, seafloor_depth):
        """Sub method of vertical_mixing, determines settling"""
        # Elements at or below seafloor are settled, by setting
        # self.elements.moving to 0.
        # These elements will not move until eventual later resuspension.
        settling = np.logical_and(self.elements.z <= seafloor_depth, self.elements.moving==1)
        if np.sum(settling) > 0:
            self.logger.debug('Settling %s elements at seafloor' % np.sum(settling))
            self.elements.moving[settling] = 0

    def resuspension(self):
        """Resuspending elements if current speed > .5 m/s"""
        resuspending = np.logical_and(self.current_speed()>.5, self.elements.moving==0)
        if np.sum(resuspending) > 0:
            # Allow moving again
            self.elements.moving[resuspending] = 1
            # Suspend 1 cm above seafloor
            self.elements.z[resuspending] = self.elements.z[resuspending] + .01
