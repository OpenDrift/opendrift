from datetime import datetime, timedelta
import numpy as np
from readers import readers

class Environment(object):

    def __init__(self):
        pass


class OpenDriftSimulation(object):

    def __init__(self, time_step=3600, proj4=None):

        # Add Readers object as single point interface
        # towards environmental data
        self.readers = readers.Readers()

        # Time step in seconds
        self.time_step = timedelta(seconds=time_step)

        if proj4 is not None:
            self.proj4 = proj4
            self.proj = readers.pyproj.Proj(self.proj4)
        else:
            self.proj4 = None

        print 'OpenDriftSimulation initialised'

    def add_reader(self, reader, **kwargs):
        self.readers.add_reader(reader, kwargs)
        if self.proj4 is None:
            self.proj4 = reader.proj4
            self.proj = readers.pyproj.Proj(self.proj4)

    def seed_point(self, lon, lat, radius, number, time, **kwargs):
        radius = radius/111000.  # convert radius from m to degrees
        kwargs['lon'] = lon + radius*(np.random.rand(number) - 0.5)
        kwargs['lat'] = lat + radius*(np.random.rand(number) - 0.5)
        self.elements = self.ElementType(**kwargs)
        self.time = time
                                          

    def run(self, steps=30):
        # - initialise / check initialisation
        # - loop:
        #   - seed particles 
        #   - call model.propagate
        #   - deactivate particles
        #   - write output
        #   - check if finished, else repeat
        self.lons = np.zeros((steps, len(self.elements)))
        self.lats = np.zeros((steps, len(self.elements)))
        for i in range(steps):
            self.get_environment()
            self.propagate()
            self.lons[i,:] = self.elements.lon
            self.lats[i,:] = self.elements.lat

    def plot(self):
        # Temporary plotting function based on Basemap
        from mpl_toolkits.basemap import Basemap
        import matplotlib.pyplot as plt
        lonmin = self.lons.min()
        lonmax = self.lons.max()
        latmin = self.lats.min()
        latmax = self.lats.max()
        buffer = 1
        map = Basemap(lonmin-buffer, latmin-buffer,
                      lonmax+buffer, latmax+buffer,
                      resolution='i', projection='cyl')
        map.plot(self.lons, self.lats, color='k')
        map.drawcoastlines()
        map.fillcontinents(color='coral')
        map.drawmeridians(np.arange(0,360,1))
        map.drawparallels(np.arange(-90,90,1))
        plt.show()

    def propagate(self):
        #print 'Getting environment data...'
        self.get_environment()
        #print 'Updating particle positions and properties...'
        self.update()
        #print 'Updating time...'
        self.time = self.time + self.time_step
        #print 'Finished'
        print self.time


    def get_environment(self):

        # Convert lon/lat to x,y
        self.elements.x, self.elements.y = self.lonlat2xy(
            self.elements.lon, self.elements.lat)
        # Find which readers have the needed parameters
        self.environment = Environment()
        environment = {}
        required_variables = self.required_variables
        for reader in self.readers.readers:
            available_variables = list(
                set(reader.variables) & set(required_variables))
            # Get available variables from this reader
            reader_x, reader_y = reader.lonlat2xy(
                self.elements.lon, self.elements.lat)  # x,y in reader coords
            env = reader.get_variables(available_variables, self.time,
                    reader_x, reader_y, self.elements.depth)
            environment.update(env)
            required_variables = list(set(required_variables) - set(available_variables))
            if len(required_variables) == 0:
                break # Got all variables, no need to check more readers

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
        self.elements.x = self.elements.x + x_vel*self.time_step.total_seconds()
        self.elements.y = self.elements.y + y_vel*self.time_step.total_seconds()
        # Calculate x,y from lon,lat
        self.elements.lon, self.elements.lat = self.xy2lonlat(
            self.elements.x, self.elements.y)

    def __repr__(self):
        outStr = '===========================\n'
        outStr += 'Model:\t' + type(self).__name__ + '\n'
        if hasattr(self, 'elements'):
            outStr += '\t%s %s particles\n' % (
                        len(self.elements), type(self.elements).__name__)
        outStr += 'Projection: %s\n' % self.proj4
        outStr += 'Readers:\n'
        for reader in self.readers.readers:
            outStr += '\t' + reader.name + '\n'
        if hasattr(self, 'time'):
            outStr += 'Time: %s\n' % (self.time)
        outStr += '===========================\n'
        return outStr
