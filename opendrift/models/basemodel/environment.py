import logging
from typing import OrderedDict, Dict, List
import copy
import traceback
import numpy as np
import pyproj

from opendrift.timer import Timeable
from opendrift.readers.basereader import BaseReader
from opendrift.readers import reader_from_url, reader_global_landmask
from opendrift.errors import NotCoveredError

from opendrift.config import Configurable

logger = logging.getLogger(__name__)


class Environment(Timeable, Configurable):
    fallback_values: Dict
    readers: OrderedDict
    priority_list: OrderedDict
    required_variables: Dict
    required_profiles_z_range: List[float]  # [min_depth, max_depth]
    max_speed: float

    proj_latlon = pyproj.Proj('+proj=latlong')

    def __init__(self, required_variables, required_profiles_z_range, max_speed, _config):
        super().__init__()

        self.fallback_values = {}
        self.readers = OrderedDict()
        self.priority_list = OrderedDict()

        self.required_variables = required_variables
        self.required_profiles_z_range = required_profiles_z_range
        self.max_speed = max_speed
        self._config = _config # reference to simulation config

        # Find variables which require profiles
        self.required_profiles = [
            var for var in self.required_variables
            if 'profiles' in self.required_variables[var]
            and self.required_variables[var]['profiles'] is True
        ]

        # Find variables which are desired, but not required
        self.desired_variables = [
            var for var in self.required_variables
            if 'important' in self.required_variables[var]
            and self.required_variables[var]['important'] is False
        ]

    def finalize(self, simulation: 'OpenDriftSimulation', simulation_extent):
        """
        Prepare environment for simulation.
        """
        self.fallback_values = simulation.get_fallback_values()
        self.__generate_constant_readers__(simulation)
        self.__add_auto_landmask__(simulation)
        self.__assert_no_missing_variables__()
        self.prepare_readers(simulation_extent, simulation.start_time,
                             simulation.expected_end_time,
                             simulation.max_speed)

    def prepare_readers(self, extent, start_time, end_time, max_speed):
        for reader in self.readers.values():
            logger.debug('\tPreparing %s' % reader.name)
            reader.prepare(extent=extent,
                           start_time=start_time,
                           end_time=end_time,
                           max_speed=max_speed)

    def __generate_constant_readers__(self, config: Configurable):
        # Make constant readers if config environment:constant:<var> is
        c = config.get_configspec('environment:constant:')
        mr = {}
        for var in list(c):
            if c[var]['value'] is not None:
                mr[var.split(':')[-1]] = c[var]['value']
        if len(mr) > 0:
            from opendrift.readers import reader_constant
            rc = reader_constant.Reader(mr)
            self.add_reader(rc, first=True)

    def __add_auto_landmask__(self, config: Configurable):
        ##############################################################
        # If no landmask has been added, we determine it dynamically
        ##############################################################
        # TODO: some more error checking here
        # If landmask is requested, it shall not be obtained from other readers
        if config.get_config('general:use_auto_landmask') is True:
            if 'land_binary_mask' in self.priority_list:
                if 'global_landmask' in self.priority_list['land_binary_mask']:
                    self.priority_list['land_binary_mask'] = [
                        'global_landmask'
                    ]
                else:
                    del self.priority_list['land_binary_mask']

        if config.get_config('general:use_auto_landmask') is True and \
                ('land_binary_mask' in self.required_variables and \
                'land_binary_mask' not in self.priority_list \
                and 'land_binary_mask' not in self.fallback_values):
            logger.info(
                'Adding a dynamical landmask with max. priority based on '
                'assumed maximum speed of %s m/s. '
                'Adding a customised landmask may be faster...' %
                self.max_speed)
            self.timer_start('preparing main loop:making dynamical landmask')
            reader_landmask = reader_global_landmask.Reader()
            self.add_reader(reader_landmask)
            self.timer_end('preparing main loop:making dynamical landmask')

    def __assert_no_missing_variables__(self):
        missing_variables = self.missing_variables()
        missing_variables = [
            m for m in missing_variables if m != 'land_binary_mask'
        ]
        if len(missing_variables) > 0:
            has_fallback = [
                var for var in missing_variables if var in self.fallback_values
            ]
            has_no_fallback = [
                var for var in missing_variables
                if var not in self.fallback_values
            ]
            #if has_fallback == missing_variables:
            if len(has_fallback) > 0:  # == missing_variables:
                logger.info('Fallback values will be used for the following '
                            'variables which have no readers: ')
                for var in has_fallback:
                    logger.info('\t%s: %f' % (var, self.fallback_values[var]))
            #else:
            if len(has_no_fallback) > 0 and len(
                    self._lazy_readers()) == 0:  # == missing_variables:
                logger.warning(
                    'No readers added for the following variables: ' +
                    str(has_no_fallback))
                raise ValueError('Readers must be added for the '
                                 'following required variables: ' +
                                 str(has_no_fallback))

    def add_readers_from_file(self, filename, timeout=10, lazy=True):
        fp = open(filename, 'r')
        sources = fp.readlines()
        sources = [line.strip() for line in sources if line[0] != '#']
        self.add_readers_from_list(sources, timeout, lazy=lazy)

    def add_readers_from_list(self,
                              urls,
                              timeout=10,
                              lazy=True,
                              variables=None):
        '''Make readers from a list of URLs or paths to netCDF datasets'''

        if isinstance(urls, str):
            urls = [urls]
        if lazy is True:
            from opendrift.readers.reader_lazy import Reader
            readers = [Reader(u) for u in urls]
            self.add_reader(readers, variables=variables)
            return

        readers = [reader_from_url(u, timeout) for u in urls]
        self.add_reader([r for r in readers if r is not None],
                        variables=variables)

    def add_reader(self, readers, variables=None, first=False):
        """Add one or more readers providing variables used by this model.

        Method may be called subsequently to add more readers
        for other variables.

        Args:
            readers: one or more (list) Reader objects.

            variables (optional): list of strings of standard_name of variables to be provided by this/these reader(s).
            first: Set to True if this reader should be set as first option
        """

        # Convert any strings to lists, for looping
        if isinstance(variables, str):
            variables = [variables]
        if isinstance(readers, BaseReader):
            readers = [readers]

        for reader in readers:
            # Check if input class is of correct type
            if not isinstance(reader, BaseReader) and \
                    not hasattr(reader, '_lazyname'):
                raise TypeError('Please provide Reader object')

            # Check that reader class contains the requested variables
            if variables is not None:
                missingVariables = set(variables) - set(reader.variables)
                if missingVariables:
                    raise ValueError(
                        'Reader %s does not provide variables: %s' %
                        (reader.name, list(missingVariables)))

            # Finally add new reader to list
            if reader.name in self.readers:
                # Reader names must be unique, adding integer
                for n in range(99999):
                    tmp_name = reader.name + '_%d' % n
                    if tmp_name not in self.readers:
                        reader.name = tmp_name
                        break

            # Horizontal buffer of reader must be large enough to cover
            # the distance possibly covered by elements within a time step
            if not reader.is_lazy:
                reader.set_buffer_size(max_speed=self.max_speed)

            self.readers[reader.name] = reader
            logger.debug('Added reader ' + reader.name)

            # Add this reader for each of the given variables
            if reader.is_lazy is False:
                for variable in variables if variables else reader.variables:
                    if variable in list(self.priority_list):
                        if reader.name not in self.priority_list[variable]:
                            if first is True:
                                self.priority_list[variable].insert(
                                    0, reader.name)
                            else:
                                self.priority_list[variable].append(
                                    reader.name)
                    else:
                        self.priority_list[variable] = [reader.name]

        # Remove/hide variables not needed by the current trajectory model
        for variable in list(self.priority_list):
            if variable not in self.required_variables:
                del self.priority_list[variable]

    def list_environment_variables(self):
        """Return list of all variables provided by the added readers."""
        variables = []
        for reader in self.readers:
            variables.extend(self.readers[reader].variables)
        return variables

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
            variables = list(self.required_variables)
        reader_groups = []
        # Find all unique reader groups
        for variable, readers in self.priority_list.items():
            if (variable in variables) and (readers not in reader_groups):
                reader_groups.append(readers)
        # Find all variables returned by the same reader group
        variable_groups = [None] * len(reader_groups)
        for variable, readers in self.priority_list.items():
            if variable not in variables:
                continue
            for i, readerGroup in enumerate(reader_groups):
                if readers == readerGroup:
                    if variable_groups[i]:
                        variable_groups[i].append(variable)
                    else:
                        variable_groups[i] = [variable]

        missing_variables = list(
            set(variables) - set(self.priority_list.keys()))

        return variable_groups, reader_groups, missing_variables

    def _lazy_readers(self):
        return [r for r in self.readers if self.readers[r].is_lazy is True]

    def _unlazy_readers(self):
        return [r for r in self.readers if self.readers[r].is_lazy is False]

    def _initialise_next_lazy_reader(self):
        '''Returns reader if successful and None if no more readers'''

        lazy_readers = self._lazy_readers()

        if len(lazy_readers) == 0:
            return None

        lazyname = lazy_readers[0]
        reader = self.readers[lazyname]

        try:
            reader.initialise()
        except Exception as e:
            logger.debug(e)
            logger.warning('Reader could not be initialised, and is'
                           ' discarded: ' + lazyname)
            self.discard_reader(reader, reason='could not be initialized')
            return self._initialise_next_lazy_reader()  # Call self

        reader.set_buffer_size(max_speed=self.max_speed)
        # Update reader lazy name with actual name
        self.readers[reader.name] = \
            self.readers.pop(lazyname)
        for var in reader.variables:
            if var in list(self.priority_list):
                self.priority_list[var].append(reader.name)
            else:
                self.priority_list[var] = [reader.name]
        # Remove variables not needed
        for variable in list(self.priority_list):
            if variable not in self.required_variables:
                del self.priority_list[variable]

        return reader

    def discard_reader_if_not_relevant(self, reader):
        if reader.is_lazy:
            return False
        if reader.start_time is not None and reader.always_valid is False:
            if hasattr(self, 'expected_end_time'
                       ) and reader.start_time > self.expected_end_time:
                self.discard_reader(reader, 'starts after simulation end')
                return True
            if hasattr(self,
                       'start_time') and reader.end_time < self.start_time:
                self.discard_reader(reader, 'ends before simuation start')
                return True
            if hasattr(self, 'time') and reader.end_time < self.time:
                self.discard_reader(reader,
                                    'ends before simuation is finished')
                return True
        if len(set(self.required_variables) & set(reader.variables)) == 0:
            self.discard_reader(
                reader, reason='does not contain any relevant variables')
            return True
        if not hasattr(reader, 'checked_for_overlap'):
            if not reader.global_coverage():
                if not hasattr(self, 'simulation_extent'):
                    logger.warning(
                        'Simulation has no simulation_extent, cannot check reader coverage'
                    )
                    return False
                # TODO
                # need a better coverage/overlap check below
                corners = reader.xy2lonlat(
                    [reader.xmin, reader.xmin, reader.xmax, reader.xmax],
                    [reader.ymax, reader.ymin, reader.ymax, reader.ymin])
                rlonmin = np.min(corners[0])
                rlonmax = np.max(corners[0])
                rlatmin = np.min(corners[1])
                rlatmax = np.max(corners[1])
                if hasattr(
                        reader, 'proj4'
                ) and 'stere' in reader.proj4 and 'lat_0=90' in reader.proj4:
                    rlatmax = 90
                if hasattr(
                        reader, 'proj4'
                ) and 'stere' in reader.proj4 and 'lat_0=-90' in reader.proj4:
                    rlatmin = -90
                if rlatmin > self.simulation_extent[3]:
                    self.discard_reader(reader, reason='too far north')
                    return True
                if rlatmax < self.simulation_extent[1]:
                    self.discard_reader(reader, reason='too far south')
                    return True
                # Disabling below checks, as +/-360 deg is not considered
                #if rlonmax < self.simulation_extent[0]:
                #    self.discard_reader(reader, reason='too far west')
                #    return True
                #if rlonmin > self.simulation_extent[2]:
                #    self.discard_reader(reader, reason='too far east')
                #    return True
            reader.checked_for_overlap = True

        return False  # Reader is not discarded

    def discard_reader(self, reader, reason):
        readername = reader.name
        logger.debug('Discarding reader (%s): %s' % (reason, readername))
        del self.readers[readername]
        if not hasattr(self, 'discarded_readers'):
            self.discarded_readers = {readername: reason}
        else:
            self.discarded_readers[readername] = reason

        # Remove from priority list
        for var in self.priority_list.copy():
            self.priority_list[var] = [
                r for r in self.priority_list[var] if r != readername
            ]
            if len(self.priority_list[var]) == 0:
                del self.priority_list[var]

    def missing_variables(self):
        """Return list of all variables for which no reader has been added."""
        return [
            var for var in self.required_variables
            if var not in self.priority_list
        ]

    def get_environment(self, variables, time, lon, lat, z, profiles):
        '''Retrieve environmental variables at requested positions.

        Updates:
            Buffer (raw data blocks) for each reader stored for performance:
                [readers].var_block_before (last before requested time)
                [readers].var_block_after (first after requested time)
                - lists of one ReaderBlock per variable group: time, x, y, [vars]

        Returns:
            environment: recarray with variables as named attributes,
                         interpolated to requested positions/time.

        '''
        self.timer_start('main loop:readers')
        print(f"{variables=}")
        # Initialise ndarray to hold environment variables
        dtype = [(var, np.float32) for var in variables]
        env = np.ma.array(np.zeros(len(lon)) * np.nan, dtype=dtype)

        if not hasattr(self, 'fallback_values'):
            self.set_fallback_values(refresh=False)

        # Discard any existing readers which are not relevant
        for readername, reader in self.readers.copy().items():
            self.discard_reader_if_not_relevant(reader)

        if 'drift:truncate_ocean_model_below_m' in self._config:
            truncate_depth = self.get_config(
                'drift:truncate_ocean_model_below_m')
            if truncate_depth is not None:
                logger.debug('Truncating ocean models below %s m' %
                             truncate_depth)
                z = z.copy()
                z[z < -truncate_depth] = -truncate_depth
                if self.required_profiles_z_range is not None:
                    self.required_profiles_z_range = np.array(
                        self.required_profiles_z_range)
                    self.required_profiles_z_range[
                        self.required_profiles_z_range <
                        -truncate_depth] = -truncate_depth

        # Initialise more lazy readers if necessary
        missing_variables = ['missingvar']
        while (len(missing_variables) > 0 and len(self._lazy_readers()) > 0):
            variable_groups, reader_groups, missing_variables = \
                self.get_reader_groups(variables)
            if hasattr(self, 'desired_variables'):
                missing_variables = list(
                    set(missing_variables) - set(self.desired_variables))
            if len(missing_variables) > 0:
                logger.debug('Variables not covered by any reader: ' +
                             str(missing_variables))
                reader = 'NotNone'
                while reader is not None:
                    reader = self._initialise_next_lazy_reader()
                    if reader is not None:
                        if self.discard_reader_if_not_relevant(reader):
                            reader = None
                    if reader is not None:
                        if (reader.covers_time(self.time) and len(
                                reader.covers_positions(lon, lat)[0]) > 0):
                            missing_variables = list(
                                set(missing_variables) - set(reader.variables))
                            if len(missing_variables) == 0:
                                break  # We cover now all variables

        # For each variable/reader group:
        variable_groups, reader_groups, missing_variables = \
            self.get_reader_groups(variables)
        for variable in variables:  # Fill with fallback value if no reader
            co = self.get_config('environment:fallback:%s' % variable)
            if co is not None:
                env[variable] = np.ma.ones(env[variable].shape) * co

        for i, variable_group in enumerate(variable_groups):
            logger.debug('----------------------------------------')
            logger.debug('Variable group %s' % (str(variable_group)))
            logger.debug('----------------------------------------')
            reader_group = reader_groups[i]
            missing_indices = np.array(range(len(lon)))
            # For each reader:
            for reader_name in reader_group:
                logger.debug('Calling reader ' + reader_name)
                logger.debug('----------------------------------------')
                self.timer_start('main loop:readers:' +
                                 reader_name.replace(':', '<colon>'))
                reader = self.readers[reader_name]
                if reader.is_lazy:
                    logger.warning('Reader is lazy, should not happen')
                    import sys
                    sys.exit('Should not happen')
                if not reader.covers_time(time):
                    logger.debug('\tOutside time coverage of reader.')
                    if reader_name == reader_group[-1]:
                        if self._initialise_next_lazy_reader() is not None:
                            logger.debug(
                                'Missing variables: calling get_environment recursively'
                            )
                            return self.get_environment(
                                variables, time, lon, lat, z, profiles)
                    continue
                # Fetch given variables at given positions from current reader
                try:
                    logger.debug('Data needed for %i elements' %
                                 len(missing_indices))
                    # Check if vertical profiles are requested from reader
                    if profiles is not None:
                        profiles_from_reader = list(
                            set(variable_group) & set(profiles))
                        if profiles_from_reader == []:
                            profiles_from_reader = None
                    else:
                        profiles_from_reader = None
                    env_tmp, env_profiles_tmp = \
                        reader.get_variables_interpolated(
                            variable_group, profiles_from_reader,
                            self.required_profiles_z_range, time,
                            lon[missing_indices], lat[missing_indices],
                            z[missing_indices], self.proj_latlon)

                except NotCoveredError as e:
                    logger.info(e)
                    self.timer_end('main loop:readers:' +
                                   reader_name.replace(':', '<colon>'))
                    if reader_name == reader_group[-1]:
                        if self._initialise_next_lazy_reader() is not None:
                            logger.debug(
                                'Missing variables: calling get_environment recursively'
                            )
                            return self.get_environment(
                                variables, time, lon, lat, z, profiles)
                    continue

                except Exception as e:  # Unknown error
                    # TODO:
                    # This could e.g. be due to corrupted files or
                    # hangig thredds-servers. A reader could be discarded
                    # after e.g. 3 such failed attempts
                    logger.info('========================')
                    logger.exception(e)
                    logger.debug(traceback.format_exc())
                    logger.info('========================')
                    self.timer_end('main loop:readers:' +
                                   reader_name.replace(':', '<colon>'))
                    if reader_name == reader_group[-1]:
                        if self._initialise_next_lazy_reader() is not None:
                            logger.debug(
                                'Missing variables: calling get_environment recursively'
                            )
                            return self.get_environment(
                                variables, time, lon, lat, z, profiles)
                    continue

                # Copy retrieved variables to env array, and mask nan-values
                for var in variable_group:
                    if var not in self.required_variables:
                        logger.debug('Not returning env-variable: ' + var)
                        continue
                    if var not in env.dtype.names:
                        continue  # Skipping variables that are only used to derive needed variables
                    env[var][missing_indices] = np.ma.masked_invalid(
                        env_tmp[var][0:len(missing_indices)]).astype('float32')
                    if profiles_from_reader is not None and var in profiles_from_reader:
                        if 'env_profiles' not in locals():
                            env_profiles = env_profiles_tmp
                        # TODO: fix to be checked
                        if var in env_profiles and var in env_profiles_tmp:
                            # If one profile has fewer vertical layers than
                            # the other, we use only the overlapping part
                            if len(env_profiles['z']) != len(
                                    env_profiles_tmp['z']):
                                logger.debug('Warning: different number of '
                                             ' vertical layers: %s and %s' %
                                             (len(env_profiles['z']),
                                              len(env_profiles_tmp['z'])))
                            z_ind = np.arange(
                                np.minimum(
                                    len(env_profiles['z']) - 1,
                                    len(env_profiles_tmp['z']) - 1))
                            # len(missing_indices) since 2 points might have been added and not removed
                            env_profiles_tmp[var] = np.ma.atleast_2d(
                                env_profiles_tmp[var])
                            env_profiles[var][np.ix_(z_ind, missing_indices)] = \
                                np.ma.masked_invalid(env_profiles_tmp[var][z_ind,0:len(missing_indices)]).astype('float32')
                            # For profiles with different numbers of layers, we extrapolate
                            if env_profiles[var].shape[0] > 1:
                                missingbottom = np.isnan(
                                    env_profiles[var][-1, :])
                                env_profiles[var][
                                    -1, missingbottom] = env_profiles[var][
                                        -2, missingbottom]

                # Detect elements with missing data, for present reader group
                if hasattr(env_tmp[variable_group[0]], 'mask'):
                    try:
                        del combined_mask
                    except:
                        pass
                    for var in variable_group:
                        tmp_var = np.ma.masked_invalid(env_tmp[var])
                        # Changed 13 Oct 2016, but uncertain of effect
                        # TODO: to be checked
                        #tmp_var = env_tmp[var]
                        if 'combined_mask' not in locals():
                            combined_mask = np.ma.getmask(tmp_var)
                        else:
                            combined_mask = \
                                np.ma.mask_or(combined_mask,
                                              np.ma.getmask(tmp_var),
                                              shrink=False)
                    try:
                        if len(missing_indices) != len(combined_mask):
                            # TODO: mask mismatch due to 2 added points
                            raise ValueError('Mismatch of masks')
                        missing_indices = missing_indices[combined_mask]
                    except Exception as ex:  # Not sure what is happening here
                        logger.info(
                            'Problems setting mask on missing_indices!')
                        logger.exception(ex)
                else:
                    missing_indices = []  # temporary workaround
                if (type(missing_indices)
                        == np.int64) or (type(missing_indices) == np.int32):
                    missing_indices = []
                self.timer_end('main loop:readers:' +
                               reader_name.replace(':', '<colon>'))
                if len(missing_indices) == 0:
                    logger.debug('Obtained data for all elements.')
                    break
                else:
                    logger.debug('Data missing for %i elements.' %
                                 (len(missing_indices)))
                    if len(self._lazy_readers()) > 0:
                        if self._initialise_next_lazy_reader() is not None:
                            logger.warning(
                                'Missing variables: calling get_environment recursively'
                            )
                            return self.get_environment(
                                variables, time, lon, lat, z, profiles)

        logger.debug('---------------------------------------')
        logger.debug('Finished processing all variable groups')

        self.timer_start('main loop:readers:postprocessing')
        for var in self.fallback_values:
            if (var not in variables) and (profiles is None
                                           or var not in profiles):
                continue
            mask = env[var].mask
            if any(mask == True):
                logger.debug(
                    '    Using fallback value %s for %s for %s elements' %
                    (self.fallback_values[var], var, np.sum(mask == True)))
                env[var][mask] = self.fallback_values[var]
            # Profiles
            if profiles is not None and var in profiles:
                if 'env_profiles' not in locals():
                    logger.debug('Creating empty dictionary for profiles not '
                                 'profided by any reader: ' +
                                 str(self.required_profiles))
                    env_profiles = {}
                    env_profiles['z'] = \
                        np.array(self.required_profiles_z_range)[::-1]
                if var not in env_profiles:
                    logger.debug(
                        '      Using fallback value %s for %s for all profiles'
                        % (self.fallback_values[var], var))
                    env_profiles[var] = self.fallback_values[var]*\
                        np.ma.ones((len(env_profiles['z']), len(lon)))
                else:
                    mask = env_profiles[var].mask
                    num_masked_values_per_element = np.sum(mask == True)
                    num_missing_profiles = np.sum(num_masked_values_per_element
                                                  == len(env_profiles['z']))
                    env_profiles[var][mask] = self.fallback_values[var]
                    logger.debug(
                        '      Using fallback value %s for %s for %s profiles'
                        % (
                            self.fallback_values[var],
                            var,
                            num_missing_profiles,
                        ))
                    num_missing_individual = np.sum(
                        num_masked_values_per_element >
                        0) - num_missing_profiles
                    if num_missing_individual > 0:
                        logger.debug(
                            '        ...plus %s individual points in other profiles'
                            % num_missing_individual)

        #######################################################
        # Some extra checks of units and realistic magnitude
        #######################################################
        if 'sea_water_temperature' in variables:
            t_kelvin = np.where(env['sea_water_temperature'] > 100)[0]
            if len(t_kelvin) > 0:
                logger.warning(
                    'Converting temperatures from Kelvin to Celcius')
                env['sea_water_temperature'][
                    t_kelvin] = env['sea_water_temperature'][t_kelvin] - 273.15
                if 'env_profiles' in locals(
                ) and 'sea_water_temperature' in env_profiles.keys():
                    env_profiles['sea_water_temperature'][:,t_kelvin] = \
                      env_profiles['sea_water_temperature'][:,t_kelvin] - 273.15

        #######################################################
        # Parameterisation of unavailable variables
        #######################################################
        if 'drift:use_tabularised_stokes_drift' in self._config and self.get_config(
                'drift:use_tabularised_stokes_drift') is True:
            if 'x_wind' not in variables:
                logger.debug('No wind available to calculate Stokes drift')
            else:
                if 'sea_surface_wave_stokes_drift_x_velocity' not in variables or (
                        env['sea_surface_wave_stokes_drift_x_velocity'].max()
                        == 0 and
                        env['sea_surface_wave_stokes_drift_y_velocity'].max()
                        == 0):
                    logger.debug('Calculating parameterised stokes drift')
                    env['sea_surface_wave_stokes_drift_x_velocity'], \
                    env['sea_surface_wave_stokes_drift_y_velocity'] = \
                        self.wave_stokes_drift_parameterised((env['x_wind'], env['y_wind']),
                            self.get_config('drift:tabularised_stokes_drift_fetch'))

                if (env['sea_surface_wave_significant_height'].max() == 0):
                    logger.debug(
                        'Calculating parameterised significant wave height')
                    env['sea_surface_wave_significant_height'] = \
                        self.wave_significant_height_parameterised((env['x_wind'], env['y_wind']),
                        self.get_config('drift:tabularised_stokes_drift_fetch'))

        #############################
        # Add uncertainty/diffusion
        #############################
        # Current
        if 'x_sea_water_velocity' in variables and \
                'y_sea_water_velocity' in variables:
            std = self.get_config('drift:current_uncertainty')
            if std > 0:
                logger.debug('Adding uncertainty for current: %s m/s' % std)
                env['x_sea_water_velocity'] += np.random.normal(
                    0, std, self.num_elements_active())
                env['y_sea_water_velocity'] += np.random.normal(
                    0, std, self.num_elements_active())
            std = self.get_config('drift:current_uncertainty_uniform')
            if std > 0:
                logger.debug('Adding uncertainty for current: %s m/s' % std)
                env['x_sea_water_velocity'] += np.random.uniform(
                    -std, std, self.num_elements_active())
                env['y_sea_water_velocity'] += np.random.uniform(
                    -std, std, self.num_elements_active())
        # Wind
        if 'x_wind' in variables and 'y_wind' in variables:
            std = self.get_config('drift:wind_uncertainty')
            if std > 0:
                logger.debug('Adding uncertainty for wind: %s m/s' % std)
                env['x_wind'] += np.random.normal(0, std,
                                                  self.num_elements_active())
                env['y_wind'] += np.random.normal(0, std,
                                                  self.num_elements_active())

        #####################
        # Diagnostic output
        #####################
        if len(env) > 0:
            logger.debug('------------ SUMMARY -------------')
            for var in variables:
                logger.debug('    %s: %g (min) %g (max)' %
                             (var, env[var].min(), env[var].max()))
            logger.debug('---------------------------------')
            # logger.debug('\t\t%s active elements' % self.num_elements_active())
            # if self.num_elements_active() > 0:
            #     lonmin = self.elements.lon.min()
            #     lonmax = self.elements.lon.max()
            #     latmin = self.elements.lat.min()
            #     latmax = self.elements.lat.max()
            #     zmin = self.elements.z.min()
            #     zmax = self.elements.z.max()
            #     if latmin == latmax:
            #         logger.debug('\t\tlatitude =  %s' % (latmin))
            #     else:
            #         logger.debug('\t\t%s <- latitude  -> %s' %
            #                      (latmin, latmax))
            #     if lonmin == lonmax:
            #         logger.debug('\t\tlongitude = %s' % (lonmin))
            #     else:
            #         logger.debug('\t\t%s <- longitude -> %s' %
            #                      (lonmin, lonmax))
            #     if zmin == zmax:
            #         logger.debug('\t\tz = %s' % (zmin))
            #     else:
            #         logger.debug('\t\t%s   <- z ->   %s' % (zmin, zmax))
            #     logger.debug('---------------------------------')

        # Prepare array indiciating which elements contain any invalid values
        missing = np.ma.masked_invalid(env[variables[0]]).mask
        for var in variables[1:]:
            missing = np.ma.mask_or(missing,
                                    np.ma.masked_invalid(env[var]).mask,
                                    shrink=False)

        # Convert dictionary to recarray and return
        if 'env_profiles' not in locals():
            env_profiles = None

        # Convert masked arrays to regular arrays for increased performance
        env = np.array(env)
        if env_profiles is not None:
            for var in env_profiles:
                env_profiles[var] = np.array(env_profiles[var])

        self.timer_end('main loop:readers:postprocessing')
        self.timer_end('main loop:readers')

        return env.view(np.recarray), env_profiles, missing

class HasEnvironment:
    """
    A class that has an `Environment`. Some shortcuts for dealing with readers are provided to the inner `env` instance.
    """
    env: Environment

    def add_reader(self, readers, variables=None, first=False):
        self.env.add_reader(readers, variables, first)
