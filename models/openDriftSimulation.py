import numpy as np
from readers import readers

class Environment(object):

    def __init__(self):
        pass


class OpenDriftSimulation(object):

    def __init__(self, time_step=None, proj4=None):

        # Add Readers object as single point interface
        # towards environmental data
        self.readers = readers.Readers()

        self.proj4 = proj4
        self.proj = readers.pyproj.Proj(self.proj4)

        print 'OpenDriftSimulation initialised'

    def seed_point(self, lon, lat, radius, number, time, **kwargs):
        radius = radius/111000.  # convert radius from m to degrees
        kwargs['lon'] = lon + radius*np.random.rand(number)
        kwargs['lat'] = lat + radius*np.random.rand(number)
        self.elements = self.ElementType(**kwargs)
        self.time = time
                                          

    def run(self):
        # - initialise / check initialisation
        # - loop:
        #   - seed particles 
        #   - call model.propagate
        #   - deactivate particles
        #   - write output
        #   - check if finished, else repeat
        self.get_environment()
        self.propagate()
        pass

    def get_environment(self):

        # Convert lon/lat to x,y
        self.elements.x, self.elements.y = self.lonlat2xy(
            self.elements.lon, self.elements.lat)
        # Find which readers have the needed parameters
        self.environment = Environment()
        environment = {}
        for reader in self.readers.readers:
            available_variables = list(
                set(reader.variables) & set(self.required_variables))
            # Get available variables from this reader
            reader_x, reader_y = reader.lonlat2xy(
                self.elements.lon, self.elements.lat)  # x,y in reader coords
            env = reader.get_variables(available_variables, self.time,
                    reader_x, reader_y, self.elements.depth)
            environment.update(env)

        for parameter in environment.keys():
            setattr(self.environment, parameter, environment[parameter])

        return self.environment

    def lonlat2xy(self, lon, lat):
        return self.proj(lon, lat, inverse=False)

    def xy2lonlat(self, x, y):
        return self.proj(x, y, inverse=True)

    def update_positions(self, x_vel, y_vel):
        # Move particles according to timestep and velocities
        # along x- and y-axes. 
        # TODO: account for vector orientation, and projection skewness

        # Calculate x,y from lon,lat
        self.elements.x, self.elements.y = self.lonlat2xy(
            self.elements.lon, self.elements.lat)
        # Update x,y
        self.elements.x = self.elements.x + x_vel*self.time_step
        self.elements.y = self.elements.y + y_vel*self.time_step
        # Calculate x,y from lon,lat
        self.elements.lon, self.elements.lat = self.xy2lonlat(
            self.elements.x, self.elements.y)

    def propagate(self):
        # Run one timestep forwards.
        # Must be overloaded by subclasses.
        pass
