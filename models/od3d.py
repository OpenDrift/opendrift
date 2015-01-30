from openDriftSimulation import OpenDriftSimulation
from elements.oil import Oil
import numpy as np

class OD3D(OpenDriftSimulation):

    ElementType = Oil
    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity']

    def update(self):
        # Use last good values, if currently missing
        try:
            mask = self.environment.x_sea_water_velocity.mask
            self.environment.x_sea_water_velocity[mask] = \
                self.environment_previous.x_sea_water_velocity[mask]
            self.environment.y_sea_water_velocity[mask] = \
                self.environment_previous.y_sea_water_velocity[mask]
        except:
            pass

        # Simply move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        # Evaporate 10% of oil mass
        evaporated = self.elements.massOil*0.1
        self.elements.massOil = self.elements.massOil - evaporated
        self.elements.massEvaporated = self.elements.massEvaporated + evaporated
