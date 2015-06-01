from openDriftSimulation import OpenDriftSimulation
from elements.oil import Oil
import numpy as np


class OpenOil(OpenDriftSimulation):
    """Open source oil trajectory model based on the OpenDrift framework.

        Developed at MET Norway based on oil weathering parameterisations
        found in open/published litterature.

        Under construction.
    """

    ElementType = Oil

    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity',
                          'sea_surface_wave_significant_height',
                          'sea_surface_wave_to_direction',
                          'x_wind', 'y_wind', 'land_binary_mask']

    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0,
                       'sea_surface_wave_significant_height': 0,
                       'sea_surface_wave_to_direction': np.nan,
                       'x_wind': 0, 'y_wind': 0}

    status_categories = {'initial': 'green',
                         'active': 'blue',
                         'dispersed': 'darkviolet',
                         'evaporated': 'yellow',
                         'stranded': 'red'}

    def update(self):
        """Update positions and properties of oil particles."""

        ## Evaporation
        windspeed = np.sqrt(self.environment.x_wind**2 +
                            self.environment.y_wind**2)

        # Evaporate 10% of oil mass
        evaporated = self.elements.massOil*0.1
        self.elements.massOil = self.elements.massOil - evaporated
        self.elements.massEvaporated = \
            self.elements.massEvaporated + evaporated

        # Dispersion
        self.elements.depth = \
            self.environment.sea_surface_wave_significant_height
        self.deactivate_elements(self.elements.depth>3.60, reason='dispersed')


        # Deactivate elements on land
        self.deactivate_elements(self.environment.land_binary_mask == 1,
                                 reason='stranded')

        # Simply move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        # Wind drag
        wind_factor = 0.02
        self.update_positions(self.environment.x_wind*wind_factor,
                              self.environment.y_wind*wind_factor)


