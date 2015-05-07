import sys
import traceback
import logging
from datetime import datetime, timedelta
from collections import OrderedDict
from abc import ABCMeta, abstractmethod, abstractproperty

import numpy as np

from readers.readers import pyproj, Reader


class OpenDriftSimulation(object):
    """Generic trajectory model class, to be extended (subclassed).

    This as an Abstract Base Class, meaning that only subclasses can
    be initiated and used.
    Any specific subclass ('model') must contain its own (or shared)
    specific type of particles (ElementType), whose properties are
    updated at each time_step using method update() on basis of model
    physics/chemistry/biology and 'required_variables' (environment)
    which are provided by one or more Reader objects.

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
        environment: recarray storing environment variables (wind, waves,
            current etc) as named attributes. Attribute names follow
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

    __metaclass__ = ABCMeta

    def __init__(self, proj4=None, seed=0, loglevel=logging.DEBUG):
        """Initialise OpenDriftSimulation

        Args:
            proj4: proj4 string defining spatial reference system.
                If not specified, SRS is taken from the first added reader.
            seed: integer or None. A given integer will yield identical
                random numbers drawn each simulation. Random numbers are
                e.g. used to distribute particles spatially when seeding,
                and may be used by modules (subclasses) for e.g. diffusion.
                Specifying a fixed value (default: 0) is useful for sensitivity
                tests. With seed = None, different random numbers will be drawn
                for subsequent runs, even with identical configuration/input.
            loglevel: set to 0 (default) to retrieve all debug information.
                Provide a higher value (e.g. 20) to receive less output.
        """
        # Dict to store readers
        self.readers = {}  # Dictionary, key=name, value=reader object
        self.priority_list = OrderedDict()
        self.use_block = True  # Set to False if reader shall perform interpolation

        # Set projection, if given
        self.set_projection(proj4)

        # Using a fixed seed will generate the same random numbers
        # each run, useful for sensitivity tests
        # Use seed = None to get different random numbers each time
        np.random.seed(seed)

        self.num_elements = 0  # Increase when seeding particles
        self.steps = 0  # Increase for each simulation step
        self.elements_deactivated = self.ElementType()  # Empty array

        logging.basicConfig(level=loglevel, format='%(levelname)s: %(message)s')

        logging.info('OpenDriftSimulation initialised')

    @abstractmethod
    def update(self):
        """Any trajectory model implementation must define an update method.
        This method must/can use environment data (self.environment) to
        update properties (including position) of its particles (self.elements)
        """

    @abstractproperty
    def ElementType(self):
        """Any trajectory model implementation must define an ElementType."""

    @abstractproperty
    def required_variables(self):
        """Any trajectory model implementation must list needed variables."""

    def set_projection(self, proj4):
        """Set the projection onto which data from readers is reprojected."""
        self.proj4 = proj4
        if proj4 is not None:
            self.proj = pyproj.Proj(self.proj4)
        else:
            self.proj = None

    def lonlat2xy(self, lon, lat):
        """Calculate x,y in own projection from given lon,lat (scalars/arrays).
        """
        if self.proj.is_latlong():
            return lon,lat
        else:
            return self.proj(lon, lat, inverse=False)

    def xy2lonlat(self, x, y):
        """Calculate lon,lat from given x,y (scalars/arrays) in own projection.
        """
        if self.proj.is_latlong():
            return x,y
        else:
            return self.proj(x, y, inverse=True)

    def add_reader(self, readers, variables=None):
        """Add one or more readers providing variables used by this model.

        Method may be called subsequently to add more readers
        for other variables.

        Args:
            readers: one or more (list) Reader objects.
            variables: optional, list of strings of standard_name of
                variables to be provided by this/these reader(s).
        """

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
                logging.debug('Added reader ' + reader.name)

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

    def missing_variables(self):
        """Return list of all variables for which no reader has been added."""
        return [var for var in self.required_variables
                   if var not in self.priority_list]

    def get_reader_groups(self, variables=None):
        """Find which groups of variables are provided by the same readers.

        This function loops through 'priority_list' (see above) and groups
        all variables returned by the same readers in the same order. This
        allows asking readers for several variables simultaneously,
        improving performance. Used by method 'get_environment'.

        Returns:
            variable_groups: list of lists of (environment) variables.
            reader_groups: list of list of reader names, corresponding to
                each of the variable_groups.
        """
        if variables is None:
            variables = self.required_variables
        reader_groups = []
        # Find all unique reader groups
        for variable, readers in self.priority_list.items():
            if readers not in reader_groups:
                reader_groups.append(readers)
        # Find all variables returned by the same reader group
        variable_groups = [None]*len(reader_groups)
        for variable, readers in self.priority_list.items():
            if variable not in variables:
                continue
            for i, readerGroup in enumerate(reader_groups):
                if readers == readerGroup:
                    if variable_groups[i]:
                        variable_groups[i].append(variable)
                    else:
                        variable_groups[i] = [variable]

        missing_variables = list(set(variables) -
                                 set(self.priority_list.keys()))

        return variable_groups, reader_groups, missing_variables

    def get_environment(self, variables, time, lon, lat, depth):
        '''Retrieve environmental variables at requested positions.

        Updates:
            Buffer (raw data blocks) for each reader stored for performace:
                [readers].var_block_before (last before requested time)
                [readers].var_block_after (first after requested time)
                    - lists of one dictionary per group:
                        - time, x, y, [vars]
        Returns:
            environment: recarray with variables as named attributes,
                         interpolated to requested positions/time.
                         Also includes "reader_number" for reference
                         to which reader is used for each element.
            
        '''
        # Initialise ndarray to hold environment variables
        dtype = [(var, np.float) for var in variables]
        env = np.ma.array(np.zeros(len(lon)), dtype=dtype)

        # For each variable/reader group:
        variable_groups, reader_groups, missing_variables = \
            self.get_reader_groups(variables)
        for variable in missing_variables:
            env[variable] = np.ma.array(env[variable].filled(
                                self.fallback_values[variable]))

        for i, variable_group in enumerate(variable_groups):
            logging.debug('Variable group %s' % (str(variable_group)))
            reader_group = reader_groups[i]
            missing_indices = np.array(range(len(lon)))
            # For each reader:
            for reader_name in reader_group:
                logging.debug('Calling reader ' + reader_name)
                reader = self.readers[reader_name]
                # Continue if not not within time
                if (reader.start_time is not None) and (
                    time < reader.start_time or time > reader.end_time):
                    logging.debug('Outside time coverage of reader (%s - %s)'
                    % (reader.start_time, reader.end_time))
                    continue
                # Fetch given variables at given positions from current reader
                try:
                    logging.debug('Data needed for %i elements' % 
                                    len(missing_indices))
                    env_tmp = reader.get_variables_from_buffer(
                        variable_group, time, lon[missing_indices], lat[missing_indices], depth, self.use_block)
                except Exception as e:
                    logging.info('========================')
                    logging.info('Exception:')
                    logging.info(e)
                    logging.debug(traceback.format_exc())
                    logging.info('========================')
                    continue
                # Rotation of vectors - TBD
                # Detect particles with missing data, continue to next reader
                    # Store reader number for the particles covered

                # Copy retrieved variables to env array
                for var in variable_group:
                    env[var][missing_indices] = env_tmp[var]

                if hasattr(env_tmp[variable_group[0]], 'mask'):
                    try:
                        del combined_mask
                    except:
                        pass
                    for var in variable_group:
                        if 'combined_mask' not in locals():
                            combined_mask = np.ma.getmask(env_tmp[var])
                        else:
                            combined_mask = np.ma.mask_or(combined_mask,
                                            np.ma.getmask(env_tmp[var]),
                                            shrink=False)
                    missing_indices = missing_indices[combined_mask]
                else:
                    missing_indices = []  # temporary workaround
                if len(missing_indices) == 0:
                    logging.debug('Obtained data for all elements.')
                    break
                else:
                    logging.debug('Data missing for %i elements' %
                                    (len(missing_indices)))
             
            # Perform default action for particles missing env data
            if len(missing_indices) > 0:
                for var in variable_group:
                    if var in self.fallback_values:
                        # Setting fallback value, presently only numeric
                        logging.debug('Using fallback value for %s: %s'
                                  % (var, self.fallback_values[var]))
                        env[var][missing_indices] = self.fallback_values[var]
        # Diagnostic output
        for var in variables:
            logging.debug('\t%s: %f (min) %f (max)' %
                            (var, env[var].min(), env[var].max()))

        # Convert dictionary to recarray and return
        return env.view(np.recarray)

    def seed_point(self, lon, lat, radius, number, time=None, **kwargs):
        """Seed a given number of particles spread around a given position.

        Arguments:
            lon: scalar, central longitude.
            lat: scalar, central latitude.
            radius: scalar, radius in meters within which particles will
                be randomly seeded.
            number: integer, number of particles to be seeded
            time: datenume, the time at which particles are seeded/released.
                If time is None (default) particles are seeded at the start
                time of the first added Reader object.
            kwargs: keyword arguments containing properties/attributes and
                values corresponding to the actual particle type (ElementType).
                These are forwarded to the ElementType class. All properties
                for which there are no default value must be specified.
        """
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
            # Use first time of first reader if time is not given for seeding
            try:
                for reader in self.readers.items():
                    if reader[1].start_time is not None:
                        firstReader = reader[1]
                        break
            except:
                raise ValueError('Time must be specified when no '
                                 'reader is added')
            logging.info('Using start time (%s) of reader %s' %
                            (firstReader.start_time, firstReader.name))
            self.time = firstReader.start_time
        else:
            self.time = time
        self.start_time = self.time  # Record start time for reference

    def deactivate_elements(self, indices):
        """Deactivate particles by moving them from elements
        to elements_deactivated."""
        # Basic, but some more housekeeping will be required later
        self.elements.move_elements(self.elements_deactivated, indices)
        logging.debug('Deactivated %i elements.' % (sum(indices)))

    def run(self, time_step=3600, steps=1000):
        """Start a trajectory simulation, after configuration.

        Performs the loop:
            - Obtain environment data for positions of all particles.
            - Call method 'update' to update (incl advect) particle properties.
        until one of the following conditions are met:
            - Maximum number of steps are reached
            - A needed variable can not be obtained by any reader
                (outside spatial/temporal domain) and has no fallback
                (default) value.
            - All particles have been deactivated (e.g. by stranding)
            - Occurance of any error, whose trace will be output to terminal.

        Before starting a model run, readers must be added for all
        reuired variables, and some particles/elements must be seeded.
        This function keeps track of the time consumed by obtaining data
        with readers, and updating their properties, respectively.

        Arguments:
            time_step: interval between particles updates, in seconds
                Default: 3600 (1 hour)
            steps: integer, maximum number of steps. End of simulation
                will be self.start_time + steps*self.time_step
        """

        missing_variables = self.missing_variables()
        if len(missing_variables) > 0:
            has_fallback = [var for var in missing_variables 
                            if var in self.fallback_values]
            if has_fallback == missing_variables:
                logging.info('Fallback values will be used for the following '
                             'variables which have no readers: ')
                for var in has_fallback:
                    logging.info('\t%s: %f' % (var, self.fallback_values[var]))
            else:
                raise ValueError('Readers must be added for the '
                                 'following required variables: '
                                 + str(missing_variables))
        # Primitive function to test overall functionality
        self.time_environment = timedelta(seconds=0)
        self.time_model = timedelta(seconds=0)
        # Time step in seconds
        self.time_step = timedelta(seconds=time_step)

        # Initialise arrays to store output (temoprary solution)
        self.lons = np.ma.array(np.zeros((steps, len(self.elements))),
                                mask=True)
        self.lats = np.ma.array(np.zeros((steps, len(self.elements))),
                                mask=True)
        for i in range(steps):
            try:
                # Get environment data
                start_time = datetime.now()
                # Display time to terminal
                logging.debug('==========================================')
                logging.info('%s - step %i of %i' % (self.time,
                             self.steps, steps))
                logging.debug('==========================================')
                self.environment = self.get_environment(
                                                self.required_variables,
                                                self.time,
                                                self.elements.lon,
                                                self.elements.lat,
                                                self.elements.depth)
                self.time_environment += datetime.now() - start_time
                start_time = datetime.now()
                # Propagate one timestep forwards
                self.update()
                self.steps += 1
                self.time_model += datetime.now() - start_time
                # Log positions
                self.lons[i, self.elements.ID-1] = self.elements.lon
                self.lats[i, self.elements.ID-1] = self.elements.lat
                if len(self.elements) == 0:
                    raise ValueError('No active elements, quitting simulation')
                logging.debug('%s active elements (%s deactivated)' %
                              (len(self.elements),
                               len(self.elements_deactivated)))
                # Updating time
                if self.time is not None:
                    self.time = self.time + self.time_step
            except Exception as e:
                logging.info('========================')
                logging.info('End of simulation:')
                logging.info(e)
                logging.info(traceback.format_exc())
                logging.info('========================')
                # Truncate lon/lat, and then return
                self.lons = self.lons[0:i-1, :]
                self.lats = self.lats[0:i-1, :]
                break

    def plot(self, background=None, buffer=.5):
        """Basic built-in plotting function intended for developing/debugging.

        Plots trajectories of all particles.
        Positions marked with colored stars:
        - green: all start positions
        - red: deactivated particles
        - blue: particles still active at end of simulation

        Requires availability of Basemap.

        Arguments:
            background: string, name of variable (standard_name) which will
                be plotted as background of trajectories, provided that it
                can be read with one of the available readers.
            buffer: float; spatial buffer of plot in degrees of
                longitude/latitude around particle collection.
        """

        if self.steps <= 2:
            logging.warning('No trajectories to plot.')
            return
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
        map.drawmeridians(np.arange(-360, 360, .5), labels=[0,0,0,0])
        map.drawparallels(np.arange(-90, 90, .5), labels=[0,1,1,0])

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
        plt.xlabel('%i active %s elements (%i deactivated)' %
                    (len(self.elements.lon), type(self.elements).__name__,
                     len(self.elements_deactivated.lon)))
        plt.title(type(self).__name__ + '  %s to %s (%i steps)' %
                    (self.start_time, self.time, self.steps))
        plt.show()

    def update_positions(self, x_vel, y_vel):
        """Move particles according to given velocity components.

        This method shall account for projection metrics (a distance
        on a map projection does not necessarily correspond to the same
        distance over true ground (not yet implemented).

        Arguments:
            x_vel and v_vel: floats, velocities in m/s of particle along
            x- and y-axes of the inherit SRS (proj4).
        """

        # Calculate x,y from lon,lat
        self.elements.x, self.elements.y = self.lonlat2xy(
            self.elements.lon, self.elements.lat)
        # Update x,y
        if self.proj.is_latlong():
            geod = pyproj.Geod(ellps='WGS84')
            velocity = np.sqrt(x_vel**2 + y_vel**2)
            azimuth = np.degrees(np.arctan2(x_vel, y_vel))  # Direction of motion
            self.elements.x, self.elements.y, back_az = geod.fwd(
                self.elements.x, self.elements.y,
                azimuth, velocity*self.time_step.total_seconds())
        else:
            self.elements.x += x_vel*self.time_step.total_seconds()
            self.elements.y += y_vel*self.time_step.total_seconds()
        # Calculate lon,lat from x,y
        self.elements.lon, self.elements.lat = self.xy2lonlat(
            self.elements.x, self.elements.y)
        # Check that new positions are valid
        if (self.elements.lon.min() < -180) or (
                self.elements.lon.min() > 360) or (
                self.elements.lat.min() < -90) or (
                self.elements.lat.max() > 90):
            logging.info('Invalid new coordinates:')
            logging.info(self.elements)
            sys.exit('Quitting')

    def __repr__(self):
        """String representation providing overview of model status."""
        outStr = '===========================\n'
        outStr += 'Model:\t' + type(self).__name__ + '\n'
        if hasattr(self, 'elements'):
            outStr += '\t%s active %s particles  (%s deactivated)\n' % (
                len(self.elements), type(self.elements).__name__,
                len(self.elements_deactivated))
        outStr += 'Projection: %s\n' % self.proj4
        variable_groups, reader_groups, missing = self.get_reader_groups()
        outStr += '-------------------\n'
        outStr += 'Environment variables:\n'
        for i, variableGroup in enumerate(variable_groups):
            outStr += '  -----\n'
            readerGroup = reader_groups[i]
            for variable in sorted(variableGroup):
                outStr += '  ' + variable + '\n'
            for i, reader in enumerate(readerGroup):
                outStr += '     ' + str(i+1) + ') ' + reader + '\n'
        if len(self.missing_variables()) > 0:
            outStr += '  -----\n'
            outStr += 'Readers not added for the following variables:\n'
            for variable in sorted(self.missing_variables()):
                outStr += '  ' + variable + '\n'
        if hasattr(self, 'time'):
            outStr += 'Time:\n'
            outStr += '\tStart: %s\n' % (self.start_time)
            outStr += '\tPresent: %s\n' % (self.time)
            if hasattr(self, 'time_step'):
                outStr += '\tSteps: %i * %s - total time: %s\n' % (
                    self.steps, self.time_step, self.time-self.start_time)
        if hasattr(self, 'time_environment'):
            outStr += 'Performance:\n'
            outStr += '\tFetching environment data: %s \n' % (
                self.time_environment)
            outStr += '\tUpdating elements: %s \n' % self.time_model
        outStr += '===========================\n'
        return outStr
