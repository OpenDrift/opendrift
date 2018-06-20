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

import logging
import numpy as np
from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements.passivetracer import PassiveTracer

# We add the property 'wind_drift_factor' to the element class
PassiveTracer.variables = PassiveTracer.add_variables([
                            ('wind_drift_factor', {'dtype': np.float32,
                                                   'unit': '%',
                                                   'default': 0.02})])


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

    def __init__(self, *args, **kwargs):

        super(OceanDrift, self).__init__(*args, **kwargs)

    #def previous_position_if(self):
    #    """Elements shall return to previous position if current is 0"""
    #    return ((self.environment.x_sea_water_velocity == 0) &
    #            (self.environment.y_sea_water_velocity == 0))

    def update(self):
        """Update positions and properties of elements."""

        # Simply move particles with ambient current
        self.advect_ocean_current()

        # Advect particles due to wind drag (according to specified
        #                                    wind_drift_factor)
        self.advect_wind()

    def seed_along_trajectory(self, lon, lat, time,
                              release_time_interval=None, **kwargs):
        '''Seed elements along given trajectory, at given time interval'''

        if release_time_interval is not None:
            # At first assuming constant time interval:
            timestep = time[1] - time[0]
            index_step = np.round(release_time_interval.total_seconds() /
                                  timestep.total_seconds())
            index_step = np.int(index_step)
            ind = np.arange(0, len(lon), index_step)
        else:
            ind = [0]  # Seed only one element at start of trajectory

        # Seed elements at given intervals
        logging.info('Seeding %s elements along trajectory' % len(ind))
        for i in ind:
            self.seed_elements(lon[i], lat[i], time=time[i], **kwargs)
