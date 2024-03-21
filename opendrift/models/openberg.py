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
# Copyright 2015, 2023, Knut-Frode Dagestad, MET Norway
# Copyright 2023, Achref Othmani, NERSC, Norway


import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray

class IcebergObj(Lagrangian3DArray):
    """Extending LagrangianArray with relevant properties for an Iceberg"""

    variables = Lagrangian3DArray.add_variables([
        ('sail', {'dtype': np.float32,	             # Sail of Iceberg (part above waterline )
                               'units': 'm',
                               'default': 10}),
        ('draft', {'dtype': np.float32,	             # Draft of Iceberg (part below waterline)
                               'units': 'm',
                               'default': 90}),
        ('length', {'dtype': np.float32,	         # length of Iceberg 
                               'units': 'm',
                               'default': 100}),
        ('width', {'dtype': np.float32,		         # width of Iceberg 
                               'units': 'm',
                               'default': 30}),
        ('weight_coeff', {'dtype': np.float32,       # Proportion of the mass of rectangular Iceberg, weight_coeff = 1.0 if 100%, 0.5 if 50% and 0.35 if 35%
                              'units': '1',
                              'default': 1}),
        ('water_drag_coeff', {'dtype': np.float32,   # cdo
                              'units': '1',
                              'default': 0.25}),
        ('wind_drag_coeff', {'dtype': np.float32,    # cda
                             'units': '1',
                             'default': 0.7})
        ])


class OpenBerg(OceanDrift):

    ElementType = IcebergObj

    # Specify which environment variables (e.g. wind, waves, currents...)
    # are needed/required by the present model, to be used for updating
    # the element properties (including propagation).
 
    required_variables = {
        'x_sea_water_velocity': {'fallback': None},
        'y_sea_water_velocity': {'fallback': None},
        'sea_surface_height': {'fallback': 0},
        'sea_floor_depth_below_sea_level': {'fallback': 10000},
        'x_wind': {'fallback': None},
        'y_wind': {'fallback': None},
        #'sea_surface_wave_significant_height': {'fallback': None},  # Needed for melting
        'land_binary_mask': {'fallback': None},
        }

    # Configuration
    def __init__(self, *args, **kwargs):

        # The constructor of parent class must always be called
        # to perform some necessary common initialisation tasks:
        super(OpenBerg, self).__init__(*args, **kwargs)


    def thermodynamics(self):
        pass  # not yet implemented
  
    def update(self):
        """Update positions and properties of icebergs"""

        self.thermodynamics()

        # Constants
        rho_water = 1027
        rho_air = 1.293
        rho_ice = 917
        rho_iceb = 900

        
        # Areas exposed     
        Ao = abs(self.elements.draft) * self.elements.length  # Area_wet
        Aa = self.elements.sail * self.elements.length        # Area_dry

        # See ACCIBERG presentation
        # https://docs.google.com/presentation/d/1O5C2v7PA3PW8a93IAGU-aS6BSKt3s-Fw/edit#slide=id.p1
        k = rho_air*self.elements.wind_drag_coeff*Aa / (rho_water*self.elements.water_drag_coeff*Ao)
        f = np.sqrt(k)/(1+np.sqrt(k))
        vx = (1-f)*self.environment.x_sea_water_velocity + f*self.environment.x_wind
        vy = (1-f)*self.environment.y_sea_water_velocity + f*self.environment.y_wind
        self.update_positions(vx, vy)

        # Grounding
        self.deactivate_elements(self.elements.draft > self.environment.sea_floor_depth_below_sea_level, reason='Grounded iceberg')
