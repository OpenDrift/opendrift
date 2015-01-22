from openDriftSimulation import OpenDriftSimulation
from elements.oil import Oil

class OD3D(OpenDriftSimulation):

    ElementType = Oil
    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity']

    def update(self):
        # Simply move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        # Evaporate 10% of oil mass
        evaporated = self.elements.massOil*0.1
        self.elements.massOil = self.elements.massOil - evaporated
        self.elements.massEvaporated = self.elements.massEvaporated + evaporated
