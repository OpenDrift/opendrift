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
# Copyright 2019, Knut-Frode Dagestad, MET Norway


import logging
import numpy as np
from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements.passivetracer import PassiveTracer

    # We add the property 'wind_drift_factor' to the element class
PassiveTracer.variables = PassiveTracer.add_variables([
                            ('wind_drift_factor', {'dtype': np.float32,
                                                'unit': '%', 
                                                'default': 0.02})])

class SeaIceDrift(OpenDriftSimulation):
    """Trajectory model based on the OpenDrift framework.

    Simply propagation with ocean sea ice (this module is not for ice bergs)
    Developed at MET Norway.

    """

    ElementType = PassiveTracer
    required_variables = [
                          'land_binary_mask', 'x_sea_water_velocity', 
                          'y_sea_water_velocity', 'x_wind', 
                          'y_wind', 'sea_ice_area_fraction',
                          'sea_ice_x_velocity', 'sea_ice_y_velocity',
                          'sea_surface_wave_significant_height',
                          'sea_surface_wave_stokes_drift_x_velocity',
                          'sea_surface_wave_stokes_drift_y_velocity']

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'x_wind': 0, 'y_wind': 0,
                       'sea_ice_area_fraction': 0,
                       'sea_ice_x_velocity' : 0,
                       'sea_ice_y_velocity' : 0,
                       'sea_surface_wave_significant_height' : 0,
                       'sea_surface_wave_stokes_drift_x_velocity' : 0,
                       'sea_surface_wave_stokes_drift_y_velocity' : 0}


    def __init__(self, *args, **kwargs):

        super(SeaIceDrift, self).__init__(*args, **kwargs)


    def advect_sea_ice(self):
        
        if hasattr(self.environment, 'sea_ice_area_fraction'):
            self.elements_on_ice(
                self.environment.sea_ice_area_fraction > 0, reason='in_ice')

        self.advect_with_sea_ice() 
        self.advect_ocean_current()  
        self.advect_wind()
        self.stokes_drift()

    def update(self):
        """Update positions and properties of elements."""

        # Move particles with sea ice velocity
        self.advect_sea_ice()