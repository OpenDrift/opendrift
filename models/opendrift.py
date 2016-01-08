# This file is part of OpenDrift.
# 
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
# 
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
# 
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import sys
import os
import glob
import types
import traceback
import logging
import pdb
from datetime import datetime, timedelta
from collections import OrderedDict
from abc import ABCMeta, abstractmethod, abstractproperty

import numpy as np
from configobj import configobj, validate
try:
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    from matplotlib import animation
except:
    logging.info('Basemap is not available, can not make plots')

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
            (presently only one implemented: seed_elements).
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
        self.elements = self.ElementType()  # Empty array

        logging.getLogger('').handlers = []
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
                    if reader.proj4 is not None and reader.proj4 != 'None':
                        self.set_projection(reader.proj4)
                        logging.debug('Using srs for common grid: %s' %
                                      self.proj4)
                    else:
                        logging.debug('%s is unprojected, cannot use '
                                      'for common grid' % reader.name)
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

        # Set projection to latlong if not taken from any of the readers
        if self.proj is None:
            self.set_projection('+proj=latlong')

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

    def get_environment(self, variables, time, lon, lat, z):
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
                        lat[missing_indices], z, self.use_block, self.proj)
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
                    logging.debug('Data missing for %i elements:' %
                                  (len(missing_indices)))

            # Perform default action for particles missing env data
            if len(missing_indices) > 0:
                for var in variable_group:
                    if var in self.fallback_values:
                        # Setting fallback value, presently only numeric
                        logging.debug('    Using fallback value for %s: %s'
                                      % (var, self.fallback_values[var]))
                        env[var][missing_indices] = self.fallback_values[var]
                    else:
                        logging.debug('\t\t%s values missing for %s' % (
                            len(missing_indices), var))
        # Diagnostic output
        logging.debug('------------ SUMMARY -------------')
        for var in variables:
            logging.debug('    %s: %g (min) %g (max)' %
                          (var, env[var].min(), env[var].max()))
        logging.debug('---------------------------------')

        # Prepare array indiciating which elements contain any invalid values
        missing = np.ma.masked_invalid(env[variables[0]]).mask
        for var in variables[1:]:
            missing = np.ma.mask_or(missing,
                                    np.ma.masked_invalid(env[var]).mask,
                                    shrink=False)

        # Convert dictionary to recarray and return
        return env.view(np.recarray), missing

    def num_elements_active(self):
        """The number of active elements."""
        if hasattr(self, 'elements'):
            return len(self.elements)
        else:
            return 0

    def num_elements_deactivated(self):
        """The number of deactivated elements."""
        if hasattr(self, 'elements_deactivated'):
            return len(self.elements_deactivated)
        else:
            return 0

    def num_elements_scheduled(self):
        if hasattr(self, 'elements_scheduled'):
            return len(self.elements_scheduled)
        else:
            return 0

    def num_elements_total(self):
        """The total number of scheduled, active and deactivated elements."""
        return self.num_elements_activated() + self.num_elements_scheduled()

    def num_elements_activated(self):
        """The total number of active and deactivated elements."""
        return self.num_elements_active() + self.num_elements_deactivated()

    def schedule_elements(self, elements, time):
        """Schedule elements to be seeded during runtime.

        Also assigns a unique ID to each particle, monotonically increasing."""

        ## TODO
        # prepare time
        if type(time) == datetime:
            time = [time]*len(elements)  # Convert to array of same length
        if not hasattr(self, 'elements_scheduled'):
            self.elements_scheduled = elements
            self.elements_scheduled_time = np.array(time)
            # We start simulation at time of release of first element:
            self.start_time = time[0]
            self.elements_scheduled.ID = np.arange(1, len(elements) + 1)
        else:
            elements.ID = np.arange(self.num_elements_scheduled() + 1,
                                    self.num_elements_scheduled() + 1 +
                                    len(elements))  # Increase ID successively
            self.elements_scheduled.extend(elements)
            self.elements_scheduled_time = np.append(
                self.elements_scheduled_time, np.array(time))


    def release_elements(self):
        """Activate elements which are scheduled within following timestep."""
        
        logging.debug('to be seeded: %s, already seeded %s' % (
                len(self.elements_scheduled), self.num_elements_activated()))
        if len(self.elements_scheduled) == 0:
            return
        if self.time_step.days >= 0: 
            indices = (self.elements_scheduled_time >= self.time) & \
                      (self.elements_scheduled_time <
                       self.time + self.time_step)
        else:
            indices = (self.elements_scheduled_time <= self.time) & \
                      (self.elements_scheduled_time >
                       self.time + self.time_step)
        self.elements_scheduled.move_elements(self.elements, indices)
        self.elements_scheduled_time = self.elements_scheduled_time[~indices]
        logging.debug('Released %i new elements.' % np.sum(indices))

    def seed_elements(self, lon, lat, radius=0, number=None, time=None,
                   cone=False, **kwargs):
        """Seed a given number of particles around given position(s).

        Arguments:
            lon: scalar or array, central longitude(s).
            lat: scalar or array, central latitude(s).
            radius: scalar or array, radius in meters around each lon-lat
                pair, within which particles will be randomly seeded.
            number: integer, total number of particles to be seeded
                Elements are spread equally among the given lon/lat points.
                Default is one particle for each lon-lat pair.
            time: datenum, the time at which particles are seeded/released.
                If time is an array with two elements, elements are seeded
                continously from start/first to end/last time.
            cone: boolean or integer. If True, lon and lat must be two element
                arrays, interpreted as the start and end position of a cone
                within which elements will be seeded. Radius may also be a 
                two element array specifying the radius around the points.
            kwargs: keyword arguments containing properties/attributes and
                values corresponding to the actual particle type (ElementType).
                These are forwarded to the ElementType class. All properties
                for which there are no default value must be specified.
        """

        if time is None:
            raise ValueError('Time of seeding must be specified')

        #################################################################
        # Make arrays of all input parameters, with one element per
        # lon/lat pair, for sequential seeding
        #################################################################
        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)
        num_points = len(lon)  # Number of lon/lat pairs
        if number is not None and number < num_points:
            raise ValueError('Number of elements must be greater or equal '
                             'to number of points.')

        if num_points > 1:
            ###############################
            # lon and lat are arrays
            ###############################
            radius_array = np.atleast_1d(radius)
            if len(radius_array) == 1:
                # If scalar radius is given, apply to all points
                radius_array = radius_array*np.ones(num_points)
            if number is None:
                # Default is one element per given position
                #number = 1*np.ones(num_points)
                number = num_points  # temporarily scalar, should be array
            number = np.atleast_1d(number)
            if len(number) == 1 and num_points > 1:
                number_array = np.floor(number/num_points)*np.ones(num_points)
                # Add remaining elements to first point
                remainder = number - np.floor(number/num_points)*num_points
                number_array[0] = number_array[0] + remainder

            if isinstance(time, datetime):
                time = [time, time]

            if type(time) == list and len(time) == 2:
                td = (time[1]-time[0])/(number-1)  # timestep between points
                time_array = [time[0] + i*td for i in range(number)]
                time_array = [t[0] for t in time_array]
                indx_time_end = np.cumsum(number_array, dtype=int)
                indx_time_start = np.append([0], indx_time_end[0:-1])
                time_array2 = [time_array[int(indx_time_start[i]):
                                          int(indx_time_end[i])]
                               for i in range(num_points)]
                time_array = time_array2  # Subset of times for this point

            if cone is True:
                ###################################################
                # lon and lat are start and end points of a cone
                ###################################################
                if len(lon) != 2 or len(lat) != 2:
                    raise ValueError('When cone is True, lon and lat must '
                                     'be 2-element arrays.')

                geod = pyproj.Geod(ellps='WGS84')
                conelonlats = geod.npts(lon[0], lat[0], lon[1], lat[1],
                                              number, radians=False)
                # Seed cone recursively
                lon, lat = zip(*conelonlats)
                if len(radius_array) == 1:
                    radius_array = [radius, radius]
                radius_array = np.linspace(radius_array[0], radius_array[1],
                                           number)
                number_array = np.ones(number)
                time_array = [time[0] + i*td for i in range(number)]

            # Recursively seeding elements around each point
            scalarargs = {}
            for i in range(len(lon)):
                for kwarg in kwargs:
                    try:
                        scalarargs[kwarg] = kwargs[kwarg][i]
                    except:
                        scalarargs[kwarg] = kwargs[kwarg]
                # Make sure to call seed function of base class,
                # not of a specific Children class
                OpenDriftSimulation.seed_elements(self,
                                lon=[lon[i]], lat=[lat[i]],
                                radius=radius_array[i],
                                number=int(number_array[i]),
                                time=time_array[i], cone=False, **scalarargs)
            return

        # Below we have only for single points
        if isinstance(time, datetime) and number > 1:
            time = [time, time]

        if type(time) == list and len(time) == 2 and number > 1:
            td = (time[1]-time[0])/(number-1)  # timestep between points
            # Spread seeding times equally between start/first and end/last
            time_array = [time[0] + i*td for i in range(number)]
        else:
            time_array = time
             
        geod = pyproj.Geod(ellps='WGS84')
        ones = np.ones(number)
        x = np.random.randn(number)*radius
        y = np.random.randn(number)*radius
        az = np.degrees(np.arctan2(x,y))
        dist = np.sqrt(x*x+y*y)
        kwargs['lon'], kwargs['lat'], az = \
            geod.fwd(lon*ones, lat*ones, az, dist, radians=False) 

        elements = self.ElementType(**kwargs)

        self.schedule_elements(elements, time_array)

    def deactivate_elements(self, indices, reason='deactivated'):
        """Schedule deactivated particles for deletion (at end of step)"""
        if sum(indices) == 0:
            return
        if reason not in self.status_categories:
            self.status_categories.append(reason)
            logging.debug('Added status %s' % (reason))
        reason_number = self.status_categories.index(reason)
        #if not hasattr(self.elements.status, "__len__"):
        if len(np.atleast_1d(self.elements.status)) == 1:
            status = np.asscalar(self.elements.status)
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
        #try:
        #    len(indices)
        #except:
        if indices == [] or len(indices) == 0:
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

        self.bufferlength = steps + 1
        # Initialise array to hold history (element properties and environment)
        # for export to file.
        #history_dtype_fields = [(name, self.elements.dtype[name]) for name in self.elements.dtype.fields]
        history_dtype_fields = [(name,
                                 self.ElementType.variables[name]['dtype'])
                                for name in self.ElementType.variables]
        # Add environment variables
        #self.history_metadata = self.elements.variables.copy()
        self.history_metadata = self.ElementType.variables.copy()
        for env_var in self.required_variables:
            history_dtype_fields.append((env_var, np.dtype('float32'))) 
            self.history_metadata[env_var] = {}
        history_dtype = np.dtype(history_dtype_fields)
        #self.history = np.ma.array(np.zeros([self.num_elements_active(),
        self.history = np.ma.array(np.zeros([len(self.elements_scheduled),
                                             self.bufferlength]),
                                   dtype=history_dtype,
                                   mask=[True])
        self.steps_exported = 0

        if outfile is not None:
            self.io_init(outfile, times=steps+1)
        else:
            self.outfile = None

        ##########################
        # Main loop
        ##########################
        for i in range(steps):
            try:
                # Get environment data
                runtime_start = datetime.now()
                # Release elements
                self.release_elements()
                # Display time to terminal
                logging.debug('==================================='*2)
                logging.info('%s - step %i of %i - %i active elements '
                             '(%i deactivated)' %
                             (self.time, self.steps + 1,
                              steps, self.num_elements_active(),
                              self.num_elements_deactivated()))
                logging.debug('%s elements scheduled.' %
                              self.num_elements_scheduled())
                logging.debug('==================================='*2)
                self.environment, missing = \
                    self.get_environment(self.required_variables,
                                         self.time,
                                         self.elements.lon,
                                         self.elements.lat,
                                         self.elements.z)

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
            logging.debug('Writing and closing output file: %s' % outfile)
            # Write buffer to outfile, and close
            if self.steps > self.steps_exported:  # Write last lines, if needed
                self.io_write_buffer()
            self.io_close()

        # Remove any elements scheduled for deactivation during last step
        self.remove_deactivated_elements()

        # Remove columns for unseeded elements in history array
        self.history = self.history[range(self.num_elements_activated()), :]
        # Remove rows for unreached timsteps in history array
        self.history = self.history[:, range(self.steps+1)]

    def state_to_buffer(self):
        """Append present state (elements and environment) to recarray."""
        # Store present state in history recarray
        for i, var in enumerate(self.elements.variables):
            # Temporarily assuming elements numbered from 0 to num_elements_active()
            # Does not hold when importing ID from a saved file, where
            # some elements have been deactivated
            self.history[var][self.elements.ID - 1,
                              self.steps - self.steps_exported] = \
                getattr(self.elements, var)
        # Copy environment data to history array
        for i, var in enumerate(self.environment.dtype.names):
            self.history[var][self.elements.ID - 1,
                              self.steps - self.steps_exported] = \
                getattr(self.environment, var)

        # Call writer if buffer is full
        if (self.outfile is not None) and \
                ((self.steps - self.steps_exported) == self.bufferlength):
            self.io_write_buffer()

    def index_of_activation_and_deactivation(self):
        """Return the indices when elements were seeded and deactivated."""

        firstlast = np.ma.notmasked_edges(self.history['lon'], axis=1)
        index_of_activation = firstlast[0][1]
        index_of_deactivation = firstlast[1][1]
        return index_of_activation, index_of_deactivation

    def set_up_map(self, buffer=.1):
        """Generate Basemap instance on which trajectories are plotted."""

        if hasattr(self, 'history'):
            lons = self.history['lon']
            lats = self.history['lat']
        else:
            lons = np.ma.array(np.reshape(self.elements.lon, (1, -1))).T
            lats = np.ma.array(np.reshape(self.elements.lat, (1, -1))).T

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

        x, y = map(lons, lats)

        try:
            firstlast = np.ma.notmasked_edges(x, axis=1)
            index_of_first = firstlast[0][1]
            index_of_last = firstlast[1][1]
        except:
            index_of_last = 0

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except:
            pass

        try:  # Maximise figure window size
            mng.resize(*mng.window.maxsize())
        except:
            pass

        return map, plt, x, y, index_of_first, index_of_last 

    def animation(self, buffer=.2, filename=None, compare=None,
                  legend=['', '']):
        """Animate last run."""

        def plot_timestep(i):
            """Sub function needed for matplotlib animation."""
            #plt.gcf().gca().set_title(str(i))
            ax.set_title(times[i])
            points.set_data(x[range(x.shape[0]), i],
                            y[range(x.shape[0]), i])
            points_deactivated.set_data(
                x_deactive[index_of_last_deactivated<i],
                y_deactive[index_of_last_deactivated<i])

            if compare is not None:
                points_other.set_data(x_other[range(x_other.shape[0]), i],
                                      y_other[range(x_other.shape[0]), i])
                points_other_deactivated.set_data(
                    x_other_deactive[index_of_last_deactivated_other<i],
                    y_other_deactive[index_of_last_deactivated_other<i])
                return points, points_other
            else:
                return points

        # Find map coordinates and plot points with empty data, to be updated
        map, plt, x, y, index_of_first, index_of_last = \
            self.set_up_map(buffer=buffer)
        ax = plt.gcf().gca()
        times = self.get_time_array()[0]
        index_of_last_deactivated = \
            index_of_last[self.elements_deactivated.ID-1]
        points = map.plot([], [], '.k', label=legend[0])[0]
        # Plot deactivated elements, with transparency
        points_deactivated = map.plot([], [], '.k', alpha=.3)[0]
        x_deactive, y_deactive = map(self.elements_deactivated.lon,
                                     self.elements_deactivated.lat)

        if compare is not None:
            if type(compare) is str:
                # Other is given as filename
                other = self.__class__(loglevel=0)
                other.io_import_file(compare)
            else:
                # Other is given as an OpenDrift object
                other = compare
            # Find map coordinates and plot data for comparison
            x_other, y_other = map(other.history['lon'], other.history['lat'])
            points_other = map.plot(x_other[0, 0], y_other[0, 0], '.r',
                                    label=legend[1])[0]
            x_other_deactive, y_other_deactive = \
                map(other.elements_deactivated.lon,
                    other.elements_deactivated.lat)
            # Plot deactivated elements, with transparency
            points_other_deactivated = map.plot([], [], '.r', alpha=.3)[0]
            firstlast = np.ma.notmasked_edges(x_other, axis=1)
            index_of_last_other = firstlast[1][1]
            index_of_last_deactivated_other = \
                index_of_last_other[other.elements_deactivated.ID-1]

        if legend != ['', '']:
            plt.legend()

        anim = animation.FuncAnimation(plt.gcf(), plot_timestep, blit=False,
                                       frames=x.shape[1], interval=50)

        if filename is not None:
            try:
                logging.info('Saving animation to ' + filename + '...')
                anim.save(filename, fps=20, clear_temp=False)
            except Exception as e:
                print 'Could not save animation:'
                logging.info(e)
                logging.debug(traceback.format_exc())

            if filename[-4:] == '.gif':
                logging.info('Making animated gif...')
                os.system('convert -delay %i _tmp*.png %s' %
                          (self.time_step.total_seconds()/3600.*24.,filename))

            logging.info('Deleting temporary figures...')
            tmp = glob.glob('_tmp*.png')
            for tfile in tmp:
                os.remove(tfile)
        else:
            plt.show()



    def plot(self, background=None, buffer=.5, linecolor=None,
             filename=None, drifter_file=None, show=True):
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

        map, plt, x, y, index_of_first, index_of_last = \
            self.set_up_map(buffer=buffer)

        # The more elements, the more transparent we make the lines
        min_alpha = 0.025
        max_elements = 5000.0
        alpha = min_alpha**(2*(self.num_elements_total()-1)/(max_elements-1))
        alpha = np.max((min_alpha, alpha))
        if hasattr(self, 'history'):
            # Plot trajectories
            if linecolor is None:
                map.plot(x.T, y.T, color='gray', alpha=alpha)
            else:
                # Color lines according to given parameter
                try:
                    param = self.history[linecolor]
                except:
                    raise ValueError('Available parameters to be used for '
                        'linecolors: ' + str(self.history.dtype.fields))
                from matplotlib.collections import LineCollection
                for i in range(x.shape[0]):
                    vind = np.arange(index_of_first[i], index_of_last[i] + 1)
                    points = np.array([x[i,vind].T, y[i,vind].T]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]],
                                              axis=1)
                    lc = LineCollection(segments,
                                        cmap=plt.get_cmap('Spectral'),
                                        norm=plt.Normalize(param.min(),
                                                           param.max()))
                    #lc.set_linewidth(3)
                    lc.set_array(param.T[vind,i])
                    plt.gca().add_collection(lc)
                axcb = plt.colorbar(lc)
                try:  # Add unit to colorbar if available
                    colorbarstring = linecolor + '  [%s]' % \
                        (self.history_metadata[linecolor]['units'])
                except:
                    colorbarstring = linecolor
                axcb.set_label(colorbarstring)

        map.scatter(x[range(x.shape[0]), index_of_first],
                    y[range(x.shape[0]), index_of_first],
                    zorder=10, edgecolor='k', linewidths=.2,
                    color=self.status_colors['initial'],
                    label='initial (%i)' % x.shape[0])
        map.scatter(x[range(x.shape[0]), index_of_last],
                    y[range(x.shape[0]), index_of_last],
                    zorder=3, edgecolor='k', linewidths=.2,
                    color=self.status_colors['active'],
                    label='active (%i)' %
                    (x.shape[0] - self.num_elements_deactivated()))

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
                           v_component[::skip, ::skip], scale=20)

        if hasattr(self, 'time'):
            plt.title(type(self).__name__ + '  %s to %s (%i steps)' %
                      (self.start_time.strftime('%Y-%m-%d %H:%M'),
                       self.time.strftime('%Y-%m-%d %H:%M'), self.steps))

        if drifter_file is not None:
            # Format of joubeh.com
            #for dfile in drifter_file:
            #    data = np.recfromcsv(dfile)
            #    x, y = map(data['longitude'], data['latitude'])
            #    map.plot(x, y, '-k', linewidth=2)
            #    map.plot(x[0], y[0], '*k')
            #    map.plot(x[-1], y[-1], '*k')

            # Format for shell buoy
            data = np.loadtxt(drifter_file, skiprows=1, usecols=(2,3))
            x, y = map(data.T[1], data.T[0])
            map.plot(x, y, '-r', linewidth=2, zorder=10)
            map.plot(x[0], y[0], '*r', zorder=10)
            map.plot(x[-1], y[-1], '*r', zorder=10)

        if filename is not None:
            #plt.savefig(filename, dpi=200)
            plt.savefig(filename)
            plt.close()
        else:
            if show is True:
                plt.show()
        
        return map, plt

    def get_time_array(self):
        """Return a list of times of last run."""
        td = self.time_step
        time_array = [self.start_time + td*i for i in range(self.steps+1)]
        time_array_relative = [td*i for i in range(self.steps+1)]
        return time_array, time_array_relative

    def plot_environment(self):
        """Plot mean wind and current velocities of element of last run."""
        x_wind = self.get_property('x_wind')[0]
        y_wind = self.get_property('x_wind')[0]
        wind = np.sqrt(x_wind**2 + y_wind**2)
        x_sea_water_velocity = self.get_property('x_sea_water_velocity')[0]
        y_sea_water_velocity = self.get_property('y_sea_water_velocity')[0]
        current = np.sqrt(x_sea_water_velocity**2 + y_sea_water_velocity**2)
        wind = np.ma.mean(wind, axis=1)
        current = np.ma.mean(current, axis=1)
        time, time_relative = self.get_time_array()
        time = np.array([t.total_seconds()/3600. for t in time_relative])

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(time, wind, 'b', label='wind speed')
        ax1.set_ylabel('Wind speed  [m/s]', color='b')
        ax1.set_xlim([0, time[-1]])
        ax1.set_ylim([0, wind.max()])

        ax2 = ax1.twinx()
        ax2.plot(time, current, 'r', label='current speed')
        ax2.set_ylabel('Current speed  [m/s]', color='r')
        ax2.set_xlim([0, time[-1]])
        ax2.set_ylim([0, current.max()])
        for tl in ax1.get_yticklabels():
            tl.set_color('b')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        ax1.set_xlabel('Time  [hours]')
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
        times = [self.start_time + n*self.time_step
                 for n in range(self.steps + 1)]
        data = self.history[prop].T[0:len(times), :]
        plt.plot(times, data)
        plt.title(prop)
        plt.xlabel('Time  [UTC]')
        try:
            plt.ylabel('%s  [%s]' %
                        (prop, self.elements.variables[prop]['units']))
        except:
            plt.ylabel(prop)
        plt.subplots_adjust(bottom=.3)
        plt.grid()
        plt.show()

    def get_property(self, propname):
        """Get property from history, sorted by status."""
        prop = self.history[propname]
        status = self.history['status']
        #for stat in self.status_categories:
        #    print '\t%s' % stat
        index_of_first, index_of_last = \
            self.index_of_activation_and_deactivation()
        j = np.arange(status.shape[1])
        # Fill arrays with last value before deactivation
        for i in range(status.shape[0]):
            status[i, j>index_of_last[i]] = status[i, index_of_last[i]]
            prop[i, j>index_of_last[i]] = prop[i, index_of_last[i]]

        return prop.T, status.T

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
        outStr += '\t%s active %s particles  (%s deactivated, %s scheduled)\n'\
            % (self.num_elements_active(), self.ElementType.__name__,
               self.num_elements_deactivated(), self.num_elements_scheduled())
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
