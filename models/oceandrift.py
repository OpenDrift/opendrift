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

from opendrift import OpenDriftSimulation
from elements.passivetracer import PassiveTracer


class OceanDrift(OpenDriftSimulation):
    """Trajectory model based on the OpenDrift framework.

    Simply propagation with ocean currents and possibly additional wind drag.
    Suitable for passive tracers or simple object like ocean drifters.
    Developed at MET Norway.

    """

    ElementType = PassiveTracer
    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity']
    required_variables.append('land_binary_mask')

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0}

    configspec = '''
        [drift]
            scheme = option('euler', 'runge-kutta', default='euler')'''

    def update(self):
        """Update positions and properties of elements."""

        # Simply move particles with ambient current
        self.advect_ocean_current()

        # Deactivate elements on land
        self.deactivate_elements(self.environment.land_binary_mask == 1,
                                 reason='stranded')
