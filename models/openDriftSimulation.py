import traceback
from datetime import datetime, timedelta
from collections import OrderedDict

import numpy as np

from readers.readers import pyproj, Reader

class Environment(object):

    def __init__(self):
        pass


class OpenDriftSimulation(object):

    def __init__(self, time_step=3600, proj4=None, seed=0):

        # Dict to store readers
        self.readers = {}  # Dictionary, key=name, value=reader object
        self.priority_list = OrderedDict()

        # Time step in seconds
        self.time_step = timedelta(seconds=time_step)

        # Set projection, if given
        self.set_projection(proj4)

        # Using a fixed seed will generate the same random numbers
        # each run, useful for sensitivity tests
        # Use seed = None to get different random numbers each time
        np.random.seed(seed)

        print 'OpenDriftSimulation initialised'

    def set_projection(self, proj4):
        self.proj4 = proj4
        if proj4 is not None:
            self.proj = pyproj.Proj(self.proj4)
        else:
            self.proj = None

    def lonlat2xy(self, lon, lat):
        return self.proj(lon, lat, inverse=False)

    def xy2lonlat(self, x, y):
        return self.proj(x, y, inverse=True)

    ######################################
    # Readers
    ######################################

    def add_reader(self, readers, variables=None):

        # Prepare lists, for looping
        if isinstance(variables, str):
            variables = [variables]
        if isinstance(readers, Reader):
            readers = [readers]

        for reader in readers:
            # Check if input class is of correct type
            if not isinstance(reader, Reader):
                raise TypeError('Please provide Reader object')

            # Check that reader class contains the requested variables
            if variables is not None:
                missingVariables = set(variables) - set(reader.variables)
                if missingVariables:
                    raise ValueError('Reader %s does not provide variables: %s' \
                        % (reader.name, list(missingVariables)))

            # Finally add new reader to list
            if reader not in self.readers:
                self.readers[reader.name] = reader
                if self.proj == None:
                    self.set_projection(reader.proj4)
                print 'Added ' + reader.name

            # Add this reader for each of the given variables
            for variable in variables if variables else reader.variables:
                if variable in self.priority_list:
                    if reader.name not in self.priority_list[variable]:
                        self.priority_list[variable].append(reader.name)
                else:
                    self.priority_list[variable] = [reader.name]
 
        # Remove/hide variables not needed by the current trajectory model
        for variable in self.priority_list:
            if variable not in self.required_variables:
                del self.priority_list[variable]

    def list_environment_variables(self):
        variables = []
        for reader in self.readers:
            variables.extend(self.readers[reader].variables)
        return variables

    def get_reader_groups(self):

        # Find which groups of variables are provided by
        # the same set of readers (in the same order)
        reader_groups = []
        # Find all unique reader groups
        for variable, readers in self.priority_list.items():
            if readers not in reader_groups:
                reader_groups.append(readers)
        # Find all variables returned by the same reader group
        variable_groups = [None]*len(reader_groups)
        for variable, readers in self.priority_list.items():
            for i, readerGroup in enumerate(reader_groups):
                if readers == readerGroup:
                    if variable_groups[i]:
                        variable_groups[i].append(variable)
                    else:
                        variable_groups[i] = [variable]
        return variable_groups, reader_groups

    ###################################################
    # Get Environment
    ###################################################
    def get_environment(self):
        '''Retrieve variables requested by this model by looping
        through all available readers'''

        if hasattr(self, 'environment'):
            self.environment_previous = self.environment  # Store last info

        # Convert lon/lat to x,y
        self.elements.x, self.elements.y = self.lonlat2xy(
            self.elements.lon, self.elements.lat)
        # Find which readers have the needed variables
        self.environment = Environment()
        environment = {}
        # Initialise environment with empty, masked arrays
        for var in self.required_variables:
            environment[var] = np.ma.array(
                np.zeros(len(self.elements.x)), mask=True)
        variable_groups, reader_groups = self.get_reader_groups()

        for i, variable_group in enumerate(variable_groups):
            reader_group = reader_groups[i]
            # Loop through readers until environment variables are found
            for reader in reader_group:
                reader = self.readers[reader]  # Use reader object, not name

                # NB: take copy, as mask ('missing') will update
                # when this element gets updated!
                missing = environment[variable_group[0]].mask.copy()
                if not True in missing:
                    print 'Found all variables, continuing...'
                    break

                # x,y in reader coords
                reader_x, reader_y = reader.lonlat2xy(
                    self.elements.lon[missing], self.elements.lat[missing])
                env = reader.get_variables(variable_group, self.time,
                        reader_x, reader_y, self.elements.depth)
                ###################################################
                # TBD: 
                # - interpolation of block onto particle array
                # - rotation of vectors
                # - addition of variable uncertainty
                ###################################################
                for var in variable_group:
                    environment[var][missing] = env[var]

        for variable in environment.keys():
            setattr(self.environment, variable, environment[variable])

        if not hasattr(self, 'environment_previous'):
            self.environment_previous = self.environment  # Use current/first

        # Use last good values, if currently missing
        #for varName in env:
        #    if varName == 'x' or varName == 'y':
        #        continue
        #    var = env[varName]
        #    mask = var.mask
        #    #self.environment.x_sea_water_velocity[mask] = \
        #    #    self.environment_previous.x_sea_water_velocity[mask]
        #    #self.environment.y_sea_water_velocity[mask] = \
        #    #    self.environment_previous.y_sea_water_velocity[mask]

    #######################
    # Run
    #######################

    def seed_point(self, lon, lat, radius, number, time, **kwargs):
        radius = radius/111000.  # convert radius from m to degrees
        kwargs['lon'] = lon + radius*(np.random.rand(number) - 0.5)
        kwargs['lat'] = lat + radius*(np.random.rand(number) - 0.5)
        self.elements = self.ElementType(**kwargs)
        if time is None:
            # Use first time of first reader of time is not given for seeding
            #print 'Using start time of reader ' + self.readers.readers[0].name
            firstReader = list(self.readers.items())[0][1]
            print 'Using start time of reader ' + firstReader.name
            self.time = firstReader.startTime
        else:
            self.time = time
        self.startTime = self.time  # Record start time for reference
                                          

    def run(self, steps=100):
        # Primitive function to test overall functionality
        self.time_environment = timedelta(seconds=0)
        self.time_model = timedelta(seconds=0)
        self.lons = np.zeros((steps, len(self.elements)))
        self.lats = np.zeros((steps, len(self.elements)))
        for i in range(steps):
            try:
                # Get environment data
                startTime = datetime.now()
                self.get_environment()
                self.time_environment += datetime.now() - startTime
                # Propagate
                startTime = datetime.now()
                self.propagate()
                self.time_model += datetime.now() - startTime
                # Log positions
                self.lons[i,:] = self.elements.lon
                self.lats[i,:] = self.elements.lat
            except Exception as e:
                print '========================'
                print 'End of simulation:'
                print e
                print traceback.format_exc()
                print '========================'
                # Truncate lon/lat, and then return
                self.lons = self.lons[0:i-1,:]
                self.lats = self.lats[0:i-1,:]
                break

    def plot(self, background=None):
        # Temporary plotting function based on Basemap
        if self.lons.shape[0] < 1:
            raise ValueError('No points to plot!')
        from mpl_toolkits.basemap import Basemap
        import matplotlib.pyplot as plt
        buffer = 1
        lonmin = self.lons.min() - buffer*2
        lonmax = self.lons.max() + buffer*2
        latmin = self.lats.min()
        latmax = self.lats.max()
        map = Basemap(lonmin-buffer, latmin-buffer,
                      lonmax+buffer, latmax+buffer,
                      resolution='h', projection='merc')
        x, y = map(self.lons, self.lats)
        map.plot(x, y, color='gray')
        map.plot(x[-1], y[-1], '*', color='r')
        map.plot(x[0], y[0], '*', color='g')
        map.drawcoastlines(color='coral')
        if background is None:
            map.fillcontinents(color='coral')
        map.drawmeridians(np.arange(0,360,1))
        map.drawparallels(np.arange(-90,90,1))
        if background is not None and None:  # Disabled
            # Plot background field, if requested
            for readerName in self.readers.readers:
                reader = self.readers.readers[readerName]
                if background in reader.variables:
                    break
                    #reader = list(self.readers.readers.items())[0][1]  # first reader
            lons, lats = map.makegrid(4, 4) # get lat/lons of ny by nx evenly space grid.
            reader_x, reader_y = reader.lonlat2xy(lons, lats)
            data = reader.get_variables(
                background, self.time-self.time_step, reader_x, reader_y,
                0, block=True)
            reader_x, reader_y = np.meshgrid(data['x'], data['y'])
            data = data[background]
            rlons, rlats = reader.xy2lonlat(reader_x, reader_y)
            map_x, map_y = map(rlons, rlats) 
            print data.shape, map_x.shape, map_y.shape
            map.contourf(map_x, map_y, data, interpolation='nearest')
        plt.show()


    def propagate(self):
        #print 'Getting environment data...'
        #self.get_environment()
        #print 'Updating particle positions and properties...'
        self.update()
        #print 'Updating time...'
        self.time = self.time + self.time_step
        #print 'Finished'
        print self.time


    def update_positions(self, x_vel, y_vel):
        # Move particles according to timestep and velocities
        # along x- and y-axes. 
        # TODO: account for vector orientation, and projection skewness

        # Calculate x,y from lon,lat
        self.elements.x, self.elements.y = self.lonlat2xy(
            self.elements.lon, self.elements.lat)
        # Update x,y
        self.elements.x += x_vel*self.time_step.total_seconds()
        self.elements.y += y_vel*self.time_step.total_seconds()
        # Calculate lon,lat from x,y
        self.elements.lon, self.elements.lat = self.xy2lonlat(
            self.elements.x, self.elements.y)

    def __repr__(self):
        outStr = '===========================\n'
        outStr += 'Model:\t' + type(self).__name__ + '\n'
        if hasattr(self, 'elements'):
            outStr += '\t%s %s particles\n' % (
                        len(self.elements), type(self.elements).__name__)
        outStr += 'Projection: %s\n' % self.proj4
        variable_groups, reader_groups = self.get_reader_groups()
        outStr += '-------------------\n'
        outStr += 'Environment variables:\n'
        for i, variableGroup in enumerate(variable_groups):
            outStr += '  -----\n'
            readerGroup = reader_groups[i]
            for variable in sorted(variableGroup):
                outStr += '  ' + variable + '\n'
            for i, reader in enumerate(readerGroup):
                outStr += '     ' + str(i+1) + ') ' + reader + '\n'
        if hasattr(self, 'time'):
            outStr += 'Time:\n'
            outStr += '\tStart: %s\n' % (self.startTime)
            outStr += '\tPresent: %s\n' % (self.time)
            outStr += '\tNumber of steps: %i\n' % (
                (self.time-self.startTime).total_seconds() /
                    self.time_step.total_seconds())
        if hasattr(self, 'time_environment'):
            outStr += 'Time spent:\n'
            outStr += '\tFetching environment data: %s \n' % self.time_environment
            outStr += '\tUpdating elements: %s \n' % self.time_model
        outStr += '===========================\n'
        return outStr
