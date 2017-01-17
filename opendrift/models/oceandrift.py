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
from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements.passivetracer import PassiveTracer

# We add the property 'wind_drift_factor' to the element class
PassiveTracer.variables = PassiveTracer.add_variables([
                            ('wind_drift_factor', {'dtype': np.float32,
                                                   'unit': '%',
                                                   'default': 0.02}),
                            ('age_seconds', {'dtype': np.float32,
                                             'units': 's',
                                             'default': 0})])


class OceanDrift(OpenDriftSimulation):
    """Trajectory model based on the OpenDrift framework.

    Simply propagation with ocean currents and possibly additional wind drag.
    Suitable for passive tracers or simple object like ocean drifters.
    Developed at MET Norway.

    """

    ElementType = PassiveTracer
    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'x_wind', 'y_wind']
    required_variables.append('land_binary_mask')

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'x_wind': 0,
                       'y_wind': 0}

    #configspec = '''
    #    [drift]
    #        scheme = option('euler', 'runge-kutta', default='euler')
    #        max_age_seconds = float(min=0, default=None)
    #        '''

    def __init__(self, *args, **kwargs):

        self._add_config('drift:max_age_seconds', 'float(min=0, default=None)',
                         'Time after which particles will be deactivated.')
        self._add_config('drift:scheme',
                         "option('euler', 'runge-kutta', default='euler')",
                         'Time after which particles will be deactivated.')

        super(OceanDrift, self).__init__(*args, **kwargs)

    def update(self):
        """Update positions and properties of elements."""

        self.elements.age_seconds += self.time_step.total_seconds()

        # Simply move particles with ambient current
        self.advect_ocean_current()

        # Advect particles due to wind drag (according to specified
        #                                    wind_drift_factor)
        self.advect_wind()

        # Deactivate elements that exceed a certain age
        if self.get_config('drift:max_age_seconds') is not None:
            self.deactivate_elements(self.elements.age_seconds >=
                                     self.get_config('drift:max_age_seconds'),
                                     reason='retired')
