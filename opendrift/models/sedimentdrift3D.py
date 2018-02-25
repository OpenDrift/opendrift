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

import numpy as np
from opendrift.models.opendrift3D import OpenDrift3DSimulation
from opendrift.models.oceandrift import OceanDrift
# from opendrift.elements.passivetracer import PassiveTracer
from opendrift.elements.buoyanttracer import BuoyantTracer
# from opendrift.elements.sedimenttracer import SedimentTracer


class SedimentDrift3D(OpenDrift3DSimulation, OceanDrift): # multiple inheritance
    """Trajectory model based on the OpenDrift framework.

    Simply propagation with horizontal and vertical ocean currents
    and possibly additional wind drag.
    Suitable for passive tracers, e.g. for tracking water particles.
    Developed at MET Norway.

    """
    # import pdb;pdb.set_trace()
    ElementType = BuoyantTracer # simply use BuoyantTracer for now - will eventually move to SedimentTracer
    required_variables = [
        'x_sea_water_velocity',
        'y_sea_water_velocity',
        'x_wind', 'y_wind',
        'upward_sea_water_velocity',
        'ocean_vertical_diffusivity',
        'sea_surface_wave_significant_height',
        'sea_surface_wave_stokes_drift_x_velocity',
        'sea_surface_wave_stokes_drift_y_velocity',
        'sea_surface_wave_period_at_variance_spectral_density_maximum',
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
        'sea_floor_depth_below_sea_level'
        ]

    required_variables.append('land_binary_mask')

    required_profiles = ['ocean_vertical_diffusivity']
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-120, 0]

    fallback_values = {
        'x_sea_water_velocity': 0,
        'y_sea_water_velocity': 0,
        'sea_surface_wave_significant_height': 0,
        'sea_surface_wave_stokes_drift_x_velocity': 0,
        'sea_surface_wave_stokes_drift_y_velocity': 0,
        'sea_surface_wave_period_at_variance_spectral_density_maximum': 0,
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 0,
        'x_wind': 0,
        'y_wind': 0,
        'upward_sea_water_velocity': 0,
        'ocean_vertical_diffusivity': 0.02,
        'sea_floor_depth_below_sea_level': 10000
        }

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(SedimentDrift3D, self).__init__(*args, **kwargs)

        # Turbulent mixing switched off by default
        self.set_config('processes:turbulentmixing', False)

    def update_terminal_velocity(self, *args, **kwargs):
        '''
        Terminal velocity due to buoyancy or sedimentation rate,
        to be used in turbulent mixing module.
        stored as an array in self.elements.terminal_velocity
        '''
        # do nothing 
        # self.elements.terminal_velocity = 0.
        pass

    def update(self):
        """Update positions and properties of elements."""

        self.elements.age_seconds += self.time_step.total_seconds()

        # Simply move particles with ambient current
        self.advect_ocean_current()

        # Advect particles due to wind drag
        # (according to specified wind_drift_factor)
        self.advect_wind()

        # Stokes drift
        self.stokes_drift()

        # Turbulent Mixing
        self.update_terminal_velocity()
        self.vertical_mixing() # using subroutine from opendrift3D.py

        # Vertical advection
        self.vertical_advection()

        # Sediment resuspension checks , if switched on
        # ToDo!

        # Deactivate elements that exceed a certain age
        if self.get_config('drift:max_age_seconds') is not None:
            self.deactivate_elements(self.elements.age_seconds >=
                                     self.get_config('drift:max_age_seconds'),
                                     reason='retired')
