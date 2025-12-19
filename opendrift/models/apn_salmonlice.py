"""
Represents salmon lice in OpenDrift.
"""

import sys
import numpy as np
from .oceandrift import OceanDrift, Lagrangian3DArray

class SalmonLouse(Lagrangian3DArray):
    """Keep track of degree days for application of mortality and larval development
    during postprocessing"""

    variables = Lagrangian3DArray.add_variables([
        ('degree_days', {'dtype': np.float32,
                         'units': 'day * degC',
                         'default': 0}),
    ]) 



class SalmonLouseDrift(OceanDrift):
    """Salmon louse trajectory model based on the OpenDrift framework.
    """

    ElementType = SalmonLouse
    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'sea_ice_area_fraction': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': 100},
        'ocean_vertical_diffusivity': {'fallback': 0.02, 'profiles': True},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
        'surface_downward_x_stress': {'fallback': 0},
        'surface_downward_y_stress': {'fallback': 0},
        'turbulent_kinetic_energy': {'fallback': 0},
        'turbulent_generic_length_scale': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
      }

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent calss
        super(SalmonLouseDrift, self).__init__(*args, **kwargs)

        # By default, eggs do not strand towards coastline
        self.set_config('general:coastline_action', 'previous')

        # Vertical mixing is enabled by default
        self.set_config('drift:vertical_mixing', True)

    def determine_vertical_velocity_A3(self):
        # Case 1 - salinity below 23 -> w = -1 mm/s
        ind1 = self.environment.sea_water_salinity < 23
        self.elements.terminal_velocity[ind1] = -1/1000

        # Case 2 - salinity between 23 and 32 -
        # linear increase in downward velocity with decreasing salinity w(S=32)=0, w(S=23)=-1
        ind2 = (self.environment.sea_water_salinity >= 23) & \
               (self.environment.sea_water_salinity <= 32)
        S = self.environment.sea_water_salinity[ind2]
        self.elements.terminal_velocity[ind2] = (1/1000) * ((1./9)*S -32./9)

        # Case 3 - salinity above 32 - w = 0.5 mm/s (upward) if daytime
        ind3 = self.environment.sea_water_salinity > 32
        if self.time.hour <= 12:
            self.elements.terminal_velocity[ind3] = 0.5/1000
        else:
            self.elements.terminal_velocity[ind3] = 0

    def update_degree_days(self):
        positive_temp_ind = self.environment.sea_water_temperature > 0
        deg_days_increase = self.environment.sea_water_temperature[positive_temp_ind] * \
                            (self.time_step.seconds/(60*60*24))
        self.elements.degree_days[positive_temp_ind] += deg_days_increase

    def update(self):

        self.determine_vertical_velocity_A3()

        self.update_degree_days()

        # Advection with horizontal ocean currents
        self.advect_ocean_current()

        # Advection with vertical ocean currents
        self.vertical_advection()

        # Terminal velocity
        self.vertical_buoyancy()
