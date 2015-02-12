from openDriftSimulation import OpenDriftSimulation
from elements.oil import Oil
import numpy as np


class OpenOil(OpenDriftSimulation):

    ElementType = Oil
    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity']
    fallvack_values = {'x_sea_water_velocity': 0,
                       'y_sea_water_velocity': 0}

    def update(self):

        # Simply move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        # Evaporate 10% of oil mass
        evaporated = self.elements.massOil*0.1
        self.elements.massOil = self.elements.massOil - evaporated
        self.elements.massEvaporated = \
            self.elements.massEvaporated + evaporated

        # Deactivate particles with no (masked) currents (on land)
        self.deactivate_elements(self.environment.x_sea_water_velocity.mask)
