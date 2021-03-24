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
ChemicalDrift is an OpenDrift module for drift and fate of chemicals.
The module is under development within the scope of the Horizon2020 project EMERGE
Manuel Aghito. Norwegian Meteorological Institute. 2020.
"""

import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.oceandrift import Lagrangian3DArray
from datetime import datetime

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

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
        'sea_surface_wave_period_at_variance_spectral_density_maximum': {'fallback': 0},
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'ocean_vertical_diffusivity': {'fallback': 0,
             'profiles': True},

        'sea_floor_depth_below_sea_level': {'fallback': 10000},
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
            logger.debug('Settling %s elements at seafloor' % np.sum(settling))
            self.elements.moving[settling] = 0

    def resuspension(self):
        """Resuspending elements if current speed > .5 m/s"""
        resuspending = np.logical_and(self.current_speed()>.5, self.elements.moving==0)
        if np.sum(resuspending) > 0:
            # Allow moving again
            self.elements.moving[resuspending] = 1
            # Suspend 1 cm above seafloor
            self.elements.z[resuspending] = self.elements.z[resuspending] + .01

    def seed_from_STEAM(self, steam, lowerbound=0, radius=0, **kwargs):
        """Seed elements based on a dataarray with STEAM emission data

        Arguments:
            steam: dataarray with steam emission data, with coordinates
                * latitude   (latitude) float32
                * longitude  (longitude) float32
                * time       (time) datetime64[ns]
            
            
            radius:      scalar, unit: meters
            lowerbound:  scalar, elements with lower values are discarded
        """

        sel=np.where(steam > lowerbound)
        t=steam.time[sel[0]].data
        la=steam.latitude[sel[1]].data
        lo=steam.longitude[sel[2]].data
        for i in range(0,t.size):
            number=np.array(steam.data[sel][i]/lowerbound).astype('int')
            self.seed_elements(lon=lo[i]*np.ones(number), lat=la[i]*np.ones(number),
                            radius=radius, number=number, time=datetime.utcfromtimestamp(t[i].astype(int) * 1e-9))
