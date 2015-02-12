from openDriftSimulation import OpenDriftSimulation
from elements.passivetracer import PassiveTracer


class WindBlow(OpenDriftSimulation):

    ElementType = PassiveTracer
    required_variables = ['x_wind', 'y_wind']
    fallvack_values = {'x_wind': 0,
                       'y_wind': 0}

    def update(self):

        # Simply move particles with ambient wind
        self.update_positions(self.environment.x_wind,
                              self.environment.y_wind)
