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

import logging; logger = logging.getLogger(__name__)
from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements.passivetracer import PassiveTracer


class WindBlow(OpenDriftSimulation):
    """Demonstration trajectory model based on OpenDrift framework.

    Simply advects a particle (passive tracer with
    no properties except for position) with the ambient wind.
    """

    ElementType = PassiveTracer
    required_variables = {
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0}
        }


    def __init__(self, *args, **kwargs):
        super(WindBlow, self).__init__(*args, **kwargs)
        self._set_config_default('drift:max_speed', 12)

    def update(self):

        # Simply move particles with ambient wind
        self.update_positions(self.environment.x_wind,
                              self.environment.y_wind)
