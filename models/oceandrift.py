from opendrift import OpenDriftSimulation
from elements.passivetracer import PassiveTracer
import numpy as np


class OceanDrift(OpenDriftSimulation):
    """Trajectory model based on the OpenDrift framework.

    Simply propagation with ocean currents.
    Developed at MET Norway for demonstration purpose.

    """

    ElementType = PassiveTracer
    required_variables = ['x_sea_water_velocity', 'y_sea_water_velocity']
    required_variables.append('land_binary_mask')

    configspec = 'runge_kutta = boolean(default=False)'

    def update(self):
        """Update positions and properties of oil particles."""

        if self.config['runge_kutta'] is True:
            self.move_rungekutta()
        else:
            # Simply move particles with ambient current
            self.update_positions(self.environment.x_sea_water_velocity,
                                  self.environment.y_sea_water_velocity)

        # Deactivate elements on land
        self.deactivate_elements(self.environment.land_binary_mask == 1,
                                 reason='stranded')

    def move_rungekutta(self):
        """Helper function for Runge-Kutta advection."""
        x_vel = self.environment.x_sea_water_velocity
        y_vel = self.environment.y_sea_water_velocity
        # Calculate x,y from lon,lat
        start_x, start_y = self.lonlat2xy(self.elements.lon, self.elements.lat)
        # Find midpoint
        mid_x = start_x + x_vel*self.time_step.total_seconds()*.5
        mid_y = start_y + y_vel*self.time_step.total_seconds()*.5
        mid_lon, mid_lat = self.xy2lonlat(mid_x, mid_y)
        # Find current at midpoint, a half timestep later
        mid_env, missing = self.get_environment(self.required_variables,
                                       self.time + self.time_step/2,
                                       mid_lon, mid_lat, self.elements.depth)

        # Move particles using runge-kutta velocity
        self.update_positions(mid_env['x_sea_water_velocity'],
                              mid_env['y_sea_water_velocity'])
