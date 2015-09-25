import sys
import types
import traceback
import logging
import pdb
from datetime import datetime, timedelta
from collections import OrderedDict
from abc import ABCMeta, abstractmethod, abstractproperty

import numpy as np
from configobj import configobj, validate

from readers.readers import pyproj, Reader


class ModelSettings(object):
    # Empty class to store model specific information,
    # to avoid namespace conflicts with OpenDriftSimulation class
    pass


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

    status_categories = ['active']  # Particles are active by default

    # Default plotting colors of trajectory endpoints
    status_colors_default = {'initial': 'green',
                             'active': 'blue',
                             'missing_data': 'gray'}

    model = ModelSettings  # To store model specific information

    configspec = ''  # Default (empty) configuration options, to be overriden

    def __init__(self, proj4=None, seed=0, iomodule='netcdf',
                 outfile=None, loglevel=logging.DEBUG):
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
            iomodule: name of module used to export data
                default: netcdf, see folder 'io' for more alternatives.
                'iomodule' is module/filename without preceeding 'io_'
            outfile: file where output from simulation is stored
            loglevel: set to 0 (default) to retrieve all debug information.
                Provide a higher value (e.g. 20) to receive less output.
        """

        # Set default configuration
        self.config = configobj.ConfigObj(configspec=self.configspec.split('\n'),
                                          raise_errors=True)
        validation = self.config.validate(validate.Validator())
        if not isinstance(validation, bool) or validation is False:
            raise ValueError('Wrong configuration: "%s"' % (validation))

        # Dict to store readers
        self.readers = {}  # Dictionary, key=name, value=reader object
        self.priority_list = OrderedDict()
        self.use_block = True  # Set to False if interpolation left to reader

        if not hasattr(self, 'fallback_values'):
            self.fallback_values = {}

        if hasattr(self, 'status_colors'):
            # Append model specific colors to (and override) default colors
            self.status_colors_default.update(self.status_colors)
            self.status_colors = self.status_colors_default
        else:
            self.status_colors = self.status_colors_default

        # Set projection, if given
        self.set_projection(proj4)

        # Using a fixed seed will generate the same random numbers
        # each run, useful for sensitivity tests
        # Use seed = None to get different random numbers each time
        np.random.seed(seed)

        self.steps = 0  # Increase for each simulation step
        self.elements_deactivated = self.ElementType()  # Empty array

        logging.basicConfig(level=loglevel,
                            format='%(levelname)s: %(message)s')

        # Prepare outfile
        try:
            io_module = __import__('export.io_' + iomodule, fromlist=
                                   ['init', 'write_buffer',
                                    'close', 'import_file'])
        except ImportError:
            logging.info('Could not import iomodule ' + iomodule)
        self.io_init = types.MethodType(io_module.init, self)
        self.io_write_buffer = types.MethodType(io_module.write_buffer, self)
        self.io_close = types.MethodType(io_module.close, self)
        self.io_import_file = types.MethodType(io_module.import_file, self)

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
            return lon, lat
        else:
            x, y = self.proj(lon, lat, inverse=False)
            if 'ob_tran' in self.proj4:
                # NB: should check if ob_tran is sufficient condition;
                # may need lonlat as well?
                return np.degrees(x), np.degrees(y)
            else:
                return x, y

    def xy2lonlat(self, x, y):
        """Calculate lon,lat from given x,y (scalars/arrays) in own projection.
        """
        if self.proj.is_latlong():
            return x, y
        else:
            if 'ob_tran' in self.proj4:
                logging.info('NB: Converting deg to rad due to ob_tran srs')
                x = np.radians(np.array(x))
                y = np.radians(np.array(y))
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
            if (variable in variables) and (readers not in reader_groups):
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
        for variable in missing_variables:  # Use fallback value, if no reader
            env[variable] = np.ma.ones(env[variable].shape)\
                * self.fallback_values[variable]

        for i, variable_group in enumerate(variable_groups):
            logging.debug('----------------------------------------')
            logging.debug('Variable group %s' % (str(variable_group)))
            logging.debug('----------------------------------------')
            reader_group = reader_groups[i]
            missing_indices = np.array(range(len(lon)))
            # For each reader:
            for reader_name in reader_group:
                logging.debug('Calling reader ' + reader_name)
                logging.debug('----------------------------------------')
                reader = self.readers[reader_name]
                # Continue if not not within time
                if (reader.start_time is not None) and (
                    (time < reader.start_time) or
                        (time > reader.end_time)):
                    logging.debug('Outside time coverage of reader '
                                  '(%s - %s)' %
                                  (reader.start_time, reader.end_time))
                    continue
                # Fetch given variables at given positions from current reader
                try:
                    logging.debug('Data needed for %i elements' %
                                  len(missing_indices))
                    env_tmp = reader.get_variables_from_buffer(
                        variable_group, time, lon[missing_indices],
                        lat[missing_indices], depth, self.use_block, self.proj)
                except Exception as e:
                    logging.info('========================')
                    logging.info('Exception:')
                    logging.info(e)
                    logging.debug(traceback.format_exc())
                    logging.info('========================')
                    continue

                # TBD: Store reader number for the particles covered

                # Copy retrieved variables to env array, and mask nan-values
                for var in variable_group:
                    env[var][missing_indices] = np.ma.masked_invalid(
                        env_tmp[var])

                # Detect elements with missing data, for present reader group
                if hasattr(env_tmp[variable_group[0]], 'mask'):
                    try:
                        del combined_mask
                    except:
                        pass
                    for var in variable_group:
                        #tmp_var = np.ma.masked_invalid(env_tmp[var])
                        tmp_var = env_tmp[var]
                        if 'combined_mask' not in locals():
                            combined_mask = np.ma.getmask(tmp_var)
                        else:
                            combined_mask = \
                                np.ma.mask_or(combined_mask,
                                              np.ma.getmask(tmp_var),
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
                    else:
                        logging.debug('%s values missing for %s' % (
                            len(missing_indices), var))
        # Diagnostic output
        for var in variables:
            logging.debug('\t%s: %f (min) %f (max)' %
                          (var, env[var].min(), env[var].max()))

        # Prepare array indiciating which elements contain any invalid values
        missing = np.ma.masked_invalid(env[variables[0]]).mask
        for var in variables[1:]:
            missing = np.ma.mask_or(missing,
                                    np.ma.masked_invalid(env[var]).mask,
                                    shrink=False)

        # Convert dictionary to recarray and return
        return env.view(np.recarray), missing

    def num_elements_active(self):
        """The number of active elements"""
        if hasattr(self, 'elements'):
            return len(self.elements)
        else:
            return 0

    def num_elements_deactivated(self):
        """The number of deactivated elements"""
        if hasattr(self, 'elements_deactivated'):
            return len(self.elements_deactivated)
        else:
            return 0

    def num_elements_total(self):
        """The total number of active and deactivated elements"""
        return self.num_elements_active() + self.num_elements_deactivated()

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
        geod = pyproj.Geod(ellps='WGS84')
        ones = np.ones(number)
        kwargs['lon'], kwargs['lat'], az = \
            geod.fwd(lon*ones, lat*ones,
                     360*np.random.rand(number),
                     radius*np.random.rand(number), radians=False)
        kwargs['ID'] = np.arange(self.num_elements_active() + 1,
                                 self.num_elements_active() + number + 1)
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
            self.start_time = firstReader.start_time
        else:
            self.start_time = time

    def deactivate_elements(self, indices, reason='deactivated'):
        """Schedule deactivated particles for deletion (at end of step)"""
        if sum(indices) == 0:
            return
        if reason not in self.status_categories:
            self.status_categories.append(reason)
            logging.debug('Added status %s' % (reason))
        reason_number = self.status_categories.index(reason)
        if not hasattr(self.elements.status, "__len__"):
            status = self.elements.status
            self.elements.status = np.zeros(self.num_elements_active())
            self.elements.status.fill(status)
        self.elements.status[indices] = reason_number
        logging.debug('%s elements scheduled for deactivation (%s)' %
                      (sum(indices), reason))

    def remove_deactivated_elements(self):
        """Moving deactivated elements from self.elements
        to self.elements_deactivated."""
        # All particles scheduled for deletion
        indices = (self.elements.status != 0)
        try:
            len(indices)
        except:
            logging.debug('No elements to deactivate')
            return  # No elements scheduled for deactivation
        # Basic, but some more housekeeping will be required later
        self.elements.move_elements(self.elements_deactivated, indices)
        logging.debug('Removed %i elements.' % (sum(indices)))
        if hasattr(self, 'environment'):
            self.environment = self.environment[~indices]
            logging.debug('Removed %i values from environment.' %
                          (sum(indices)))
            #if self.num_elements_active() == 0:
            #    raise ValueError('No more active elements.')  # End simulation

    def run(self, time_step=3600, steps=1000, outfile=None):
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
        required variables, and some particles/elements must be seeded.
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

        # Some cleanup needed if starting from imported state
        if self.steps > 1:
            self.steps = 0
        if hasattr(self, 'history'):
            # Delete history matrix before new run
            delattr(self, 'history')
            # Renumbering elements from 0 to num_elements, necessary fix when
            # importing from file, where elements may have been deactivated
            self.elements.ID = np.arange(0, self.num_elements_active())

        # Store runtime to report on OpenDrift performance
        self.runtime_environment = timedelta(seconds=0)
        self.runtime_model = timedelta(seconds=0)

        self.time_step = timedelta(seconds=time_step)  # Time step in seconds
        self.time = self.start_time  # Start time has been set when seeding

        if outfile is not None:
            self.io_init(outfile, times=steps+1)
        else:
            self.outfile = None
        self.bufferlength = steps + 1

        for i in range(steps):
            try:
                # Get environment data
                runtime_start = datetime.now()
                # Display time to terminal
                logging.debug('==================================='*2)
                logging.info('%s - step %i of %i - %i active elements '
                             '(%i deactivated)' %
                             (self.time, self.steps + 1,
                              steps, self.num_elements_active(),
                              self.num_elements_deactivated()))
                logging.debug('==================================='*2)
                self.environment, missing = \
                    self.get_environment(self.required_variables,
                                         self.time,
                                         self.elements.lon,
                                         self.elements.lat,
                                         self.elements.depth)

                self.deactivate_elements(missing, reason='missing_data')

                self.runtime_environment += datetime.now() - runtime_start
                runtime_start = datetime.now()
                self.state_to_buffer()  # Append status to outfile
                #if self.steps > 0:
                self.remove_deactivated_elements()
                # Propagate one timestep forwards
                self.steps += 1
                if self.num_elements_active() == 0:
                    raise ValueError('No more active elements, quitting.')
                logging.debug('Calling module: %s.update()' %
                              type(self).__name__)
                self.update()
                self.runtime_model += datetime.now() - runtime_start
                if self.num_elements_active() == 0:
                    raise ValueError('No active elements, quitting simulation')
                logging.debug('%s active elements (%s deactivated)' %
                              (self.num_elements_active(),
                               self.num_elements_deactivated()))
                # Updating time
                if self.time is not None:
                    self.time = self.time + self.time_step
            except Exception as e:
                logging.info('========================')
                logging.info('End of simulation:')
                logging.info(e)
                logging.info(traceback.format_exc())
                logging.info('========================')
                break

        logging.debug('Cleaning up')

        self.state_to_buffer()  # Append final status to buffer

        if outfile is not None:
            # Write buffer to outfile, and close
            if self.steps > self.steps_exported:  # Write last lines, if needed
                self.io_write_buffer()
            self.io_close()

        # Remove any elements scheduled for deactivation during last step
        self.remove_deactivated_elements()

    def state_to_buffer(self):
        """Append present state (elements and environment) to recarray."""
        if not hasattr(self, 'history'):
            # Presently only element properties are stored, env TBD
            self.history = np.ma.array(np.zeros([self.num_elements_active(),
                                       self.bufferlength]),
                                       dtype=self.elements.dtype,
                                       mask=[True])
            self.steps_exported = 0
        # Store present state in history recarray
        for i, var in enumerate(self.elements.variables):
            # Temporarily assuming elements numbered from 0 to num_elements_active()
            # Does not hold when importing ID from a saved file, where
            # some elements have been deactivated
            self.history[var][self.elements.ID - 1,
                              self.steps - self.steps_exported] = \
                getattr(self.elements, var)

        # Call writer if buffer is full
        if (self.outfile is not None) and \
                ((self.steps - self.steps_exported) == self.bufferlength):
            self.io_write_buffer()

    def plot(self, background=None, buffer=.5,
             filename=None, drifter_file=None):
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

        if hasattr(self, 'history'):
            lons = self.history['lon']
            lats = self.history['lat']
        else:
            lons = np.ma.array(np.reshape(self.elements.lon, (1, -1))).T
            lats = np.ma.array(np.reshape(self.elements.lat, (1, -1))).T

        try:
            from mpl_toolkits.basemap import Basemap
            import matplotlib.pyplot as plt
        except:
            sys.exit('Basemap is needed to plot trajectories')

        # Initialise map
        lonmin = lons.min() - buffer*2
        lonmax = lons.max() + buffer*2
        latmin = lats.min() - buffer
        latmax = lats.max() + buffer
        if 'basemap_landmask' in self.readers:
            # Using an eventual Basemap already used to check stranding
            map = self.readers['basemap_landmask'].map
        else:
            # Otherwise create a new Basemap covering the elements
            map = Basemap(lonmin, latmin, lonmax, latmax,
                          resolution='h', projection='merc')

        map.drawcoastlines(color='gray')
        if background is None:  # Fill continents if no field is requested
            map.fillcontinents(color='#ddaa99')
        # Adjusting spacing of lon-lat lines dynamically
        latspan = map.latmax - map.latmin
        if latspan > 20:
            deltalat = 1
        elif latspan > 1 and latspan < 20:
            deltalat = .5
        else:
            deltalat = .1
        map.drawmeridians(np.arange(np.floor(map.lonmin),
                                    np.ceil(map.lonmax), deltalat),
                          labels=[0, 0, 0, 1])
        map.drawparallels(np.arange(np.floor(map.latmin),
                                    np.ceil(map.latmax), deltalat),
                          labels=[0, 1, 1, 0])

        # Trajectories
        x, y = map(lons, lats)
        # The more elements, the more transparent we make the lines
        min_alpha = 0.025
        max_elements = 5000.0
        alpha = min_alpha**(2*(self.num_elements_total()-1)/(max_elements-1))
        alpha = np.max((min_alpha, alpha))
        if hasattr(self, 'history'):
            map.plot(x.T, y.T, color='gray', alpha=alpha)  # Plot trajectories
        # NB: should check that endpoints are always included

        map.scatter(x[:, 0], y[:, 0], zorder=10, edgecolor='k', linewidths=.2,
                    color=self.status_colors['initial'],
                    label='initial (%i)' % x.shape[0])
        map.scatter(x[:, -1], y[:, -1], zorder=3, edgecolor='k', linewidths=.2,
                    color=self.status_colors['active'],
                    label='active (%i)' %
                    (x.shape[0] - self.num_elements_deactivated()))

        try:
            index_of_last = (~x.mask).sum(axis=1) - 1
        except:
            index_of_last = 0
        x_deactivated, y_deactivated = map(self.elements_deactivated.lon,
                                           self.elements_deactivated.lat)
        # Plot deactivated elements, labeled by deactivation reason
        for statusnum, status in enumerate(self.status_categories):
            if status == 'active':
                continue  # plotted above
            if status not in self.status_colors:
                # If no color specified, pick an unused one
                for color in ['red', 'blue', 'green', 'black', 'gray',
                              'cyan', 'DarkSeaGreen', 'brown']:
                    if color not in self.status_colors.values():
                        self.status_colors[status] = color
                        break
            indices = np.where(self.elements_deactivated.status == statusnum)
            if len(indices[0]) > 0:
                map.scatter(x_deactivated[indices], y_deactivated[indices],
                            zorder=3, edgecolor='k', linewidths=.1,
                            color=self.status_colors[status],
                            label='%s (%i)' % (status, len(indices[0])))
        try:
            plt.legend()
        except Exception as e:
            print 'Cannot plot legend, due to bug in matplotlib:'
            print traceback.format_exc()

        if background is not None:  # Disabled
            # Plot background field, if requested
            if type(background) is list:
                variable = background[0]  # A vector is requested
            else:
                variable = background  # A scalar is requested
            for readerName in self.readers:
                reader = self.readers[readerName]
                if variable in reader.variables:
                    break
            # Get lat/lons of ny by nx evenly space grid.
            lons, lats = map.makegrid(4, 4)
            reader_x, reader_y = reader.lonlat2xy(lons, lats)
            data = reader.get_variables(
                background, self.time-self.time_step, reader_x, reader_y,
                0, block=True)
            reader_x, reader_y = np.meshgrid(data['x'], data['y'])
            if type(background) is list:
                u_component = data[background[0]]
                v_component = data[background[1]]
                scalar = np.sqrt(u_component**2 + v_component**2)
                # NB: rotation not completed!
                u_component, v_component = reader.rotate_vectors(
                    reader_x, reader_y, u_component, v_component,
                    reader.proj4, map.srs)
            else:
                scalar = data[background]
            rlons, rlats = reader.xy2lonlat(reader_x, reader_y)
            map_x, map_y = map(rlons, rlats)
            map.pcolormesh(map_x, map_y, scalar, alpha=1)

            if type(background) is list:
                skip = 10
                map.quiver(map_x[::skip, ::skip], map_y[::skip, ::skip],
                           u_component[::skip, ::skip],
                           v_component[::skip, ::skip], scale=10)

        if hasattr(self, 'time'):
            plt.title(type(self).__name__ + '  %s to %s (%i steps)' %
                      (self.start_time.strftime('%Y-%m-%d %H:%M'),
                       self.time.strftime('%Y-%m-%d %H:%M'), self.steps))

        if drifter_file is not None:
            for dfile in drifter_file:
                data = np.recfromcsv(dfile)
                x, y = map(data['longitude'], data['latitude'])
                map.plot(x, y, '-k', linewidth=2)
                map.plot(x[0], y[0], '*k')
                map.plot(x[-1], y[-1], '*k')

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except:
            pass

        try:  # Maximise figure window size
            mng.resize(*mng.window.maxsize())
        except:
            pass

        if filename is not None:
            #plt.savefig(filename, dpi=200)
            plt.savefig(filename)
            plt.close()
        else:
            plt.show()

    def plot_property(self, prop):
        """Basic function to plot time series of any element properties."""
        import matplotlib.pyplot as plt
        from matplotlib import dates

        hfmt = dates.DateFormatter('%d %b %Y %H:%M')
        fig = plt.figure()
        ax = fig.gca()
        ax.xaxis.set_major_formatter(hfmt)
        plt.xticks(rotation='vertical')
        data = self.history[prop].T
        times = [self.start_time + n*self.time_step
                 for n in range(self.steps + 1)]
        plt.plot(times, data)
        plt.title(prop)
        plt.xlabel('Time  [UTC]')
        plt.ylabel('%s  [%s]' % (prop, self.elements.variables[prop]['units']))
        plt.subplots_adjust(bottom=.3)
        plt.grid()
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

        geod = pyproj.Geod(ellps='WGS84')

        azimuth = np.degrees(np.arctan2(x_vel, y_vel))  # Direction of motion
        velocity = np.sqrt(x_vel**2 + y_vel**2)  # Velocity in m/s

        if not self.proj.is_latlong():  # Need to rotate SRS
            # Calculate x,y from lon,lat
            self.elements.x, self.elements.y = self.lonlat2xy(
                self.elements.lon, self.elements.lat)
            # Calculate azimuth orientation of y-axis at particle locations
            delta_y = 1000  # Using delta of 1000 m to calculate azimuth
            lon2, lat2 = self.xy2lonlat(self.elements.x,
                                        self.elements.y + delta_y)
            azimuth_srs = geod.inv(self.elements.lon, self.elements.lat,
                                   lon2, lat2)[0]
            azimuth = azimuth + azimuth_srs

        # Calculate new positions
        self.elements.lon, self.elements.lat, back_az = geod.fwd(
            self.elements.lon, self.elements.lat,
            azimuth, velocity*self.time_step.total_seconds())

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
        outStr += '\t%s active %s particles  (%s deactivated)\n' % (
            self.num_elements_active(), self.ElementType.__name__,
            self.num_elements_deactivated())
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
        if hasattr(self, 'runtime_environment'):
            outStr += 'Performance:\n'
            outStr += '\tFetching environment data: %s \n' % (
                self.runtime_environment)
            outStr += '\tUpdating elements: %s \n' % self.runtime_model
        outStr += '===========================\n'
        return outStr
