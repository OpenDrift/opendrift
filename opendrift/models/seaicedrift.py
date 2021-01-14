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
# Copyright 2019, Knut-Frode Dagestad, MET Norway

import logging; logger = logging.getLogger(__name__)
from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements.passivetracer import PassiveTracer


class SeaIceDrift(OpenDriftSimulation):
    """Trajectory model based on the OpenDrift framework.

    Simply propagation with ocean sea ice (this module is not for ice bergs)
    Developed at MET Norway.

    """

    ElementType = PassiveTracer
    required_variables = {
            'sea_ice_x_velocity': {'fallback': None},
            'sea_ice_y_velocity': {'fallback': None},
            'land_binary_mask': {'fallback': None}
        }

    def __init__(self, *args, **kwargs):

        super(SeaIceDrift, self).__init__(*args, **kwargs)


    def update(self):
        """Update positions and properties of elements."""

        # Move particles with sea ice velocity
        self.advect_with_sea_ice()
