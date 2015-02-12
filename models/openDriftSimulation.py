import sys
import traceback
from datetime import datetime, timedelta
from collections import OrderedDict

import numpy as np

from readers.readers import pyproj, Reader


class Environment(object):
    """Store environment variables (from readers) as named attributes."""

    def __init__(self):
        pass


class OpenDriftSimulation(object):
    """Generic trajectory model class, to be extended (subclassed) by modules.

    Attributes:
        ElementType: the type (class) of particles to be used by this model
        elements: object of the class ElementType, storing the specific
            particle properties (ndarrays and scalars) of all active particles
            as named attributes. Elements are added by seeding-functions
            (presently only one implemented: seed_point).
        elements_deactivated: ElementType object containing particles which
            have been deactivated (and removed from 'elements')
        required_variables: list of strings of CF standard_names which is
            needed by this model (update function) to update properties of
            particles ('elements') at each time_step. This core class has
            no required_elements, this is implemented by subclasses/modules.
        environment: Environment object storing environment variables (wind,
            waves, current etc) as named attributes. Attribute names follow
            standard_name from CF-convention, allowing any OpenDriftSimulation
            module/subclass using environment data from any readers which
            can provide the requested variables. Used in method 'update'
            to update properties of elements every time_step.
        proj4: string defining the common spatial reference system (SRS) onto
            which data from all readers are interpolated
        proj: Proj object initialised from proj4 string; used for
            coordinate tranformations
        time_step: timedelta object, time interval at which element properties
            are updated (including advection).
        readers: Dictionary where values are Reader objects, and names are
            unique reference keywords used to access a given reader (typically
            filename or URL)
        priority_list: OrderedDict where names are variable names,
            and values are lists of names (kewywords) of the reader, in the
            order of priority (user defined) of which the readers shall be
            called to retrieve environmental data.
    """

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

        self.num_elements = 0  # Increase when seeding particles
        self.elements_deactivated = self.ElementType()  # Empty array

        print 'OpenDriftSimulation initialised'

    def set_projection(self, proj4):
        """Set the projection onto which data from readers is reprojected."""
        self.proj4 = proj4
        if proj4 is not None:
            self.proj = pyproj.Proj(self.proj4)
        else:
            self.proj = None

    def lonlat2xy(self, lon, lat):
        return self.proj(lon, lat, inverse=False)

    def xy2lonlat(self, x, y):
        return self.proj(x, y, inverse=True)

    def add_reader(self, readers, variables=None):
        """Add one or more readers providing variables used by this model."""

        # Convert any strings to lists, for looping
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
                    raise ValueError('Reader %s does not provide variables: %s'
                                     % (reader.name, list(missingVariables)))

            # Finally add new reader to list
            if reader not in self.readers:
                self.readers[reader.name] = reader
                if self.proj is None:
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
        """Return list of all variables provided by the added readers."""
        variables = []
        for reader in self.readers:
            variables.extend(self.readers[reader].variables)
        return variables

    def get_reader_groups(self):
        """Find which groups of variables are provided by the same readers.
        """
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
                    break

                # x,y may be different for each reader, related to reader.proj4
                reader_x, reader_y = reader.lonlat2xy(
                    self.elements.lon[missing], self.elements.lat[missing])
                env = reader.get_variables(variable_group, self.time,
                                           reader_x, reader_y,
                                           self.elements.depth)
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

    #######################
    # Run
    #######################

    def seed_point(self, lon, lat, radius, number, time=None, **kwargs):
        radius = radius/111000.  # convert radius from m to degrees
        kwargs['lon'] = lon + radius*(np.random.rand(number) - 0.5)
        kwargs['lat'] = lat + radius*(np.random.rand(number) - 0.5)
        kwargs['ID'] = np.arange(self.num_elements + 1,
                                 self.num_elements + number + 1)
        if hasattr(self, 'elements'):
            self.elements.extend(self.ElementType(**kwargs))
        else:
            self.elements = self.ElementType(**kwargs)
        if time is None:
            # Use first time of first reader of time is not given for seeding
            try:
                firstReader = list(self.readers.items())[0][1]
            except:
                raise ValueError('Time must be specified when no reader is added')
            print 'Using start time of reader ' + firstReader.name
            self.time = firstReader.startTime
        else:
            self.time = time
        self.startTime = self.time  # Record start time for reference

    def deactivate_elements(self, indices):
        ''' Basic, but some other housekeeping may be needed later '''
        self.elements.move_elements(self.elements_deactivated, indices)

    def run(self, steps=1000000):
        # Primitive function to test overall functionality
        self.time_environment = timedelta(seconds=0)
        self.time_model = timedelta(seconds=0)
        # Initialise arrays to store output (temoprary solution)
        self.lons = np.ma.array(np.zeros((steps, len(self.elements))),
                                mask=True)
        self.lats = np.ma.array(np.zeros((steps, len(self.elements))),
                                mask=True)
        for i in range(steps):
            try:
                # Get environment data
                startTime = datetime.now()
                self.get_environment()
                self.time_environment += datetime.now() - startTime
                startTime = datetime.now()
                # Propagate one timestep forwards
                self.update()
                # Updating time
                self.time = self.time + self.time_step
                # Display time to terminal
                print self.time
                self.time_model += datetime.now() - startTime
                # Log positions
                self.lons[i, self.elements.ID-1] = self.elements.lon
                self.lats[i, self.elements.ID-1] = self.elements.lat
                if len(self.elements) == 0:
                    raise ValueError('No active elements, quitting simulation')
            except Exception as e:
                print '========================'
                print 'End of simulation:'
                print e
                print traceback.format_exc()
                print '========================'
                # Truncate lon/lat, and then return
                self.lons = self.lons[0:i-1, :]
                self.lats = self.lats[0:i-1, :]
                break

    def plot(self, background=None, buffer=.5):
        # Temporary plotting function based on Basemap
        if self.lons.shape[0] < 1:
            raise ValueError('No points to plot!')
        try:
            from mpl_toolkits.basemap import Basemap
            import matplotlib.pyplot as plt
        except:
            sys.exit('Basemap is needed to plot trajectories')

        # Initialise map
        lonmin = self.lons.min() - buffer*2
        lonmax = self.lons.max() + buffer*2
        latmin = self.lats.min()
        latmax = self.lats.max()
        map = Basemap(lonmin-buffer, latmin-buffer,
                      lonmax+buffer, latmax+buffer,
                      resolution='h', projection='merc')
        map.drawcoastlines(color='gray')
        if background is None:  # Fill continents if no field is requested
            map.fillcontinents(color='#ddaa99')
        map.drawmeridians(np.arange(0, 360, .5))
        map.drawparallels(np.arange(-90, 90, .5))

        # Trajectories
        x, y = map(self.lons, self.lats)
        numParticles = x.shape[1]
        index_of_last = (~x.mask).sum(axis=0)-1
        map.plot(x, y, color='gray')  # Plot trajectories
        map.plot(x[index_of_last, range(numParticles)],  # Deactivated
                 y[index_of_last, range(numParticles)], '*', color='r')
        map.plot(x[-1], y[-1], '*', color='b')  # Active
        map.plot(x[0], y[0], '*', color='g')  # Seed positions
        if background is not None:  # Disabled
            # Plot background field, if requested
            for readerName in self.readers:
                reader = self.readers[readerName]
                if background in reader.variables:
                    break
            # Get lat/lons of ny by nx evenly space grid.
            lons, lats = map.makegrid(4, 4)
            reader_x, reader_y = reader.lonlat2xy(lons, lats)
            data = reader.get_variables(
                background, self.time-self.time_step, reader_x, reader_y,
                0, block=True)
            reader_x, reader_y = np.meshgrid(data['x'], data['y'])
            data = data[background]
            rlons, rlats = reader.xy2lonlat(reader_x, reader_y)
            map_x, map_y = map(rlons, rlats)
            map.contourf(map_x, map_y, data, interpolation='nearest')
        plt.show()

    def propagate(self):
        # Updating particle positions and properties
        self.update()
        # Updating time
        self.time = self.time + self.time_step
        # Display time
        print self.time

    def update_positions(self, x_vel, y_vel):
        # Move particles according to timestep and velocities
        # along x- and y-axes.
        # TODO: account for vector orientation and projection metrics

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
            outStr += '\t%s active %s particles  (%s deactivated)\n' % (
                len(self.elements), type(self.elements).__name__,
                len(self.elements_deactivated))
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
            outStr += '\tFetching environment data: %s \n' % (
                self.time_environment)
            outStr += '\tUpdating elements: %s \n' % self.time_model
        outStr += '===========================\n'
        return outStr
