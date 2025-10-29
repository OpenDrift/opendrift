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
# Copyright 2022, Knut-Frode Dagestad, MET Norway

"""
OpenHNS is a 3D HNS drift module bundled within the OpenDrift framework.
"""

import numpy as np
import logging

logger = logging.getLogger(__name__)

from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED


# Defining the oil element properties
class HNS(Lagrangian3DArray):
    """Extending LagrangianArray with variables relevant for HNS particles."""

    variables = Lagrangian3DArray.add_variables([
        (   'mass', {
                'dtype': np.float32,
                'units': 'kg',
                'seed': False,
                'default': 1
        }),
         (   'mass_evaporated', {
                'dtype': np.float32,
                'units': 'kg',
                'seed': False,
                'default': 0
        }),
         (   'mass_dissolved', {
                'dtype': np.float32,
                'units': 'kg',
                'seed': False,
                'default': 0
        }),
        (
            'viscosity',
            {
                'dtype': np.float32,
                'units': 'N s/m2 (Pa s)',
                'seed': False,
                'default': 0.005
            }),
        (
            'density',
            {
                'dtype': np.float32,
                'units': 'kg/m^3',
                'seed': False,
                'default': 880
            }),
        (
            'wind_drift_factor',
            {
                'dtype':
                np.float32,  # TODO: inherit from
                'units':
                '%',  # OceanDrift
                'description':
                'Elements at the ocean surface are moved by '
                'this fraction of the wind vector, in addition to '
                'currents and Stokes drift',
                'default':
                0.02
            }),
        (
           'diameter',
            {
                'dtype': np.float32,  # Droplet diameter
                'units': 'm',
                'seed': False,
                'default': 0.
            })
    ])


class OpenHNS(OceanDrift):
    """Open source HNS drift model based on the OpenDrift framework.

        Developed at MET Norway based on parameterisations
        found in open/published litterature.

        Under construction.
    """

    ElementType = HNS

    required_variables = {
        'x_sea_water_velocity': {
            'fallback': None
        },
        'y_sea_water_velocity': {
            'fallback': None
        },
        'sea_surface_height': {'fallback': 0},
        'x_wind': {
            'fallback': None
        },
        'y_wind': {
            'fallback': None
        },
        'upward_sea_water_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_significant_height': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_stokes_drift_x_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_stokes_drift_y_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_period_at_variance_spectral_density_maximum': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment':
        {
            'fallback': 0,
            'important': False
        },
        'sea_ice_area_fraction': {
            'fallback': 0,
            'important': False
        },
        'sea_ice_x_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_ice_y_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_water_temperature': {
            'fallback': 10,
            'profiles': True
        },
        'sea_water_salinity': {
            'fallback': 34,
            'profiles': True
        },
        'sea_floor_depth_below_sea_level': {
            'fallback': 10000
        },
        'ocean_vertical_diffusivity': {
            'fallback': 0.02,
            'important': False,
            'profiles': True
        },
        'land_binary_mask': {
            'fallback': None
        },
        'ocean_mixed_layer_thickness': {
            'fallback': 50,
            'important': False
        },
    }


    hns_types = {
        'butyl': {'evaporation_rate': .03, 'dissolution_rate': .05},
        'acetone': {'evaporation_rate': .16, 'dissolution_rate': .01},
        'xylene': {'evaporation_rate': .25, 'dissolution_rate': .1}
        }

    # Default colors for plotting
    status_colors = {
        'initial': 'green',
        'active': 'blue',
        'missing_data': 'gray',
        'stranded': 'red',
        'evaporated': 'yellow',
        'dispersed': 'magenta'
    }


    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(OpenHNS, self).__init__(*args, **kwargs)

        self._add_config({
           'seed:hns_type': {
                'type':
                'enum',
                'enum':
                list(self.hns_types),
                'default':
                list(self.hns_types)[0],
                'level':
                CONFIG_LEVEL_ESSENTIAL,
                'description':
                'HNS type to be used for the simulation'
            },
        })

        self._set_config_default('drift:vertical_advection', False)
        self._set_config_default('drift:vertical_mixing', True)
        self._set_config_default('drift:current_uncertainty', 0.05)
        self._set_config_default('drift:wind_uncertainty', 0.5)

    def seed_elements(self, hns_type=None, *args, **kwargs):
        if hns_type is not None:
            self.set_config('seed:hns_type', hns_type)
        self.hns_type = self.hns_types[self.get_config('seed:hns_type')]
        super(OpenHNS, self).seed_elements(*args, **kwargs)

    def evaporation(self):
        surface = np.where(self.elements.z == 0)[0]
        random_number = np.random.uniform(0, 1, len(surface))
        evaporated = np.where(random_number > 1-self.hns_type['evaporation_rate'])[0]
        logger.debug('Evaporating %i of %i elements at ocean surface' % (len(evaporated), len(surface)))
        self.elements.wind_drift_factor[surface[evaporated]] = 1  # Shall follow wind 100%
        self.elements.z[surface[evaporated]] = 10  # Moving evaporated elements to 10m height
        self.elements.mass_evaporated[surface[evaporated]] = self.elements.mass[surface[evaporated]]
        self.elements.mass[surface[evaporated]] = 0

    def dissolution(self):
        surface = np.where(self.elements.z == 0)[0]
        random_number = np.random.uniform(0, 1, len(surface))
        dissolved = np.where(random_number > 1-self.hns_type['dissolution_rate'])[0]
        logger.debug('Dissolving %i of %i elements at ocean surface' % (len(dissolved), len(surface)))
        self.elements.wind_drift_factor[surface[dissolved]] = 0  # Submerged, no windage
        self.elements.z[surface[dissolved]] = -10  # Moving dissolved elements to 10m depth
        self.elements.mass_dissolved[surface[dissolved]] = self.elements.mass[surface[dissolved]]
        self.elements.mass[surface[dissolved]] = 0

    def update(self):

        self.evaporation()
        self.dissolution()
        self.update_terminal_velocity()
        #self.vertical_mixing()
        self.advect_ocean_current()
        self.stokes_drift()
        self.advect_wind()
