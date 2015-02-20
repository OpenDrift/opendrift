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
                          'x_wind', 'y_wind']
    fallback_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0}

    def update(self):
        """Update positions and properties of oil particles."""

        # Simply move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        # Wind drag
        wind_factor = 0.02
        self.update_positions(self.environment.x_wind*wind_factor,
                              self.environment.y_wind*wind_factor)

        ## Dispersion
        #windspeed = np.sqrt(self.environment.x_wind**2 +
        #                    self.environment.y_wind**2)
        #wave_height = 0.0246*windspeed**2  # Martinsen et al. (1994)
        #at_surface = np.where(self.elements.depth == 0)

        # Evaporate 10% of oil mass
        evaporated = self.elements.massOil*0.1
        self.elements.massOil = self.elements.massOil - evaporated
        self.elements.massEvaporated = \
            self.elements.massEvaporated + evaporated

        # Deactivate particles with no (masked) currents (on land)
        #self.deactivate_elements(self.environment.x_sea_water_velocity.mask)
