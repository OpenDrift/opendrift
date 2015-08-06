from opendrift import OpenDriftSimulation
from elements.passivetracer import PassiveTracer


class WindBlow(OpenDriftSimulation):
    """Demonstration trajectory model based on OpenDrift framework.

    Simply advects a particle (passive tracer with
    no properties except for position) with the ambient wind.
    """

    ElementType = PassiveTracer
    required_variables = ['x_wind', 'y_wind']
    fallback_values = {'x_wind': 0,
                       'y_wind': 0}

    def update(self):

        # Simply move particles with ambient wind
        self.update_positions(self.environment.x_wind,
                              self.environment.y_wind)
