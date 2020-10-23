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
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway
# Copyright 2020, Gaute Hope, MET Norway

import sys
import logging
import copy
from bisect import bisect_left
from abc import abstractmethod, ABCMeta
from datetime import datetime, timedelta
from collections import OrderedDict
from multiprocessing import Process, Manager, cpu_count

from scipy.interpolate import LinearNDInterpolator
from scipy.ndimage import map_coordinates
import numpy as np

from opendrift.readers.interpolation import ReaderBlock

import pyproj

try:
    basestring
except NameError:
    basestring = str

# Som valid (but extreme) ranges for checking that values are reasonable
standard_names = {
    'x_wind': {'valid_min': -50, 'valid_max': 50},
    'y_wind': {'valid_min': -50, 'valid_max': 50},
    'x_sea_water_velocity': {'valid_min': -10, 'valid_max': 10},
    'y_sea_water_velocity': {'valid_min': -10, 'valid_max': 10},
    'ocean_vertical_diffusivity': {'valid_min': 0, 'valid_max': 1}}

# Identify x-y vector components/pairs for rotation (NB: not east-west pairs!)
vector_pairs_xy = [
    ['x_wind', 'y_wind'],
    ['sea_ice_x_velocity', 'sea_ice_y_velocity'],
    ['x_sea_water_velocity', 'y_sea_water_velocity'],
    ['sea_surface_wave_stokes_drift_x_velocity',
     'sea_surface_wave_stokes_drift_y_velocity']
    ]


class fakeproj():
    # For readers with unprojected domain, we emulate a
    # pyproj class with needed functions
    class _crs:
        is_geographic = False

    crs = _crs()

    def __call__(self, x, y, inverse=False):
        # Simply return x and y since these are also row/column indices
        return x, y

    srs = 'None'


class BaseReader(object):
    """Parent Reader class, to be subclassed by specific readers.
    """

    __metaclass__ = ABCMeta

    return_block = True  # By default, all readers should be
                         # capable of returning blocks of data

    # Default interpolation method, see function interpolate_block()
    interpolation = 'linearNDFast'

    verticalbuffer = 1  # To be overridden by application as needed

    start_time = None

    # Mapping variable names, e.g. from east-north to x-y, temporarily
    # presuming coordinate system then is lon-lat for equivalence
    variable_aliases = {
        'sea_water_potential_temperature': 'sea_water_temperature',
        'x_wind_10m': 'x_wind',
        'y_wind_10m': 'y_wind',
        'sea_water_x_velocity': 'x_sea_water_velocity',
        'sea_water_y_velocity': 'y_sea_water_velocity',
        'x_sea_ice_velocity': 'sea_ice_x_velocity',
        'y_sea_ice_velocity': 'sea_ice_y_velocity',
        'barotropic_sea_water_x_velocity': 'sea_ice_x_velocity',
        'barotropic_sea_water_y_velocity': 'sea_ice_y_velocity',
        'salinity_vertical_diffusion_coefficient' : 'ocean_vertical_diffusivity'
        }

    xy2eastnorth_mapping = {
        'x_sea_water_velocity': ['eastward_sea_water_velocity', 'surface_eastward_sea_water_velocity',
                                 'eastward_current_velocity', 'eastward_tidal_current',
                                 'eastward_ekman_current_velocity', 'eastward_geostrophic_current_velocity',
                                 'eastward_eulerian_current_velocity', 'surface_geostrophic_eastward_sea_water_velocity',
                                 'surface_geostrophic_eastward_sea_water_velocity_assuming_sea_level_for_geoid',
                                 'surface_eastward_geostrophic_sea_water_velocity_assuming_sea_level_for_geoid'],
        'y_sea_water_velocity': ['northward_sea_water_velocity', 'surface_northward_sea_water_velocity',
                                 'northward_current_velocity', 'northward_tidal_current',
                                 'northward_ekman_current_velocity', 'northward_geostrophic_current_velocity',
                                 'northward_eulerian_current_velocity', 'surface_geostrophic_northward_sea_water_velocity',
                                 'surface_geostrophic_northward_sea_water_velocity_assuming_sea_level_for_geoid',
                                 'surface_northward_geostrophic_sea_water_velocity_assuming_sea_level_for_geoid'],
        'x_wind': 'eastward_wind', 'y_wind': 'northward_wind'}

    logger = logging.getLogger('opendrift')  # using common logger

    def __init__(self):
        # Common constructor for all readers

        # Dictionaries to store blocks of data for reuse (buffering)
        self.var_block_before = {}  # Data for last timestep before present
        self.var_block_after = {}   # Data for first timestep after present

        self.always_valid = False  # Set to True if a single field should
                                   # be valid at all times

        self.is_lazy = False  # Generally False

        # Set projection for coordinate transformations
        self.simulation_SRS = False  # Avoid unnecessary vector rotation
        if hasattr(self, 'proj'):
            self.projected = True
        else:
            if hasattr(self, 'proj4'):
                self.projected = True
                try:
                    self.proj = pyproj.Proj(self.proj4)
                except:
                    # Workaround for proj-issue with zero flattening:
                    # https://github.com/OSGeo/proj.4/issues/1191
                    origproj4 = self.proj4
                    self.proj4 = self.proj4.replace('+e=0.0', '')
                    self.proj4 = self.proj4.replace('+e=0', '')
                    self.proj4 = self.proj4.replace('+f=0.0', '')
                    self.proj4 = self.proj4.replace('+f=0', '')
                    if origproj4 != self.proj4:
                        self.logger.info('Removing flattening parameter from proj4; %s -> %s' % (origproj4, self.proj4))
                    self.proj = pyproj.Proj(self.proj4)
            else:
                self.proj4 = 'None'
                self.proj = fakeproj()
                self.projected = False
                self.logger.info('Making Splines for lon,lat to x,y conversion...')
                self.xmin = self.ymin = 0.
                self.delta_x = self.delta_y = 1.
                self.xmax = self.lon.shape[1] - 1
                self.ymax = self.lon.shape[0] - 1
                self.numx = self.xmax
                self.numy = self.ymax
                block_x, block_y = np.meshgrid(
                    np.arange(self.xmin, self.xmax + 1, 1),
                    np.arange(self.ymin, self.ymax + 1, 1))

                # Making interpolator (lon, lat) -> x
                self.spl_x = LinearNDInterpolator((self.lon.ravel(),
                                                   self.lat.ravel()),
                                                  block_x.ravel(),
                                                  fill_value=np.nan)
                # Reusing x-interpolator (deepcopy) with data for y
                self.spl_y = copy.deepcopy(self.spl_x)
                self.spl_y.values[:, 0] = block_y.ravel()
                # Call interpolator to avoid threading-problem:
                # https://github.com/scipy/scipy/issues/8856
                self.spl_x((0,0)), self.spl_y((0,0))

        # Check if there are holes in time domain
        if not hasattr(self, 'time_step'):
            self.time_step = None
        if not hasattr(self, 'start_time'):
            self.start_time = None
        if not hasattr(self, 'end_time'):
            self.end_time = None
        if self.start_time is not None and self.time_step is not None:# and len(self.times) > 1:
            self.expected_time_steps = (
                self.end_time - self.start_time).total_seconds() / (
                self.time_step.total_seconds()) + 1
            if hasattr(self, 'times'):
                self.missing_time_steps = self.expected_time_steps - \
                    len(self.times)
            else:
                self.missing_time_steps = 0
            self.actual_time_steps = self.expected_time_steps - \
                self.missing_time_steps

        # Making sure start_time is datetime, and not cftime object
        if self.start_time is not None:
             self.start_time = datetime(self.start_time.year, self.start_time.month,
                                   self.start_time.day, self.start_time.hour,
                                   self.start_time.minute, self.start_time.second)

        # Calculate shape (size) of domain
        try:
            numx = (self.xmax - self.xmin)/self.delta_x + 1
            numy = (self.ymax - self.ymin)/self.delta_y + 1
            self.shape = (int(numx), int(numy))
        except:
            self.shape = None

        self.set_buffer_size(max_speed=5)  # To be overriden by user/model

        # Check if there are east/north-oriented vectors
        for var in self.variables:
            for xvar, eastnorthvar in self.xy2eastnorth_mapping.items():
                if xvar in self.variables:
                    continue  # We have both x/y and east/north components
                if var in eastnorthvar:
                    self.logger.info('Variable %s will be rotated from %s' % (xvar, var))
                    self.variables.append(xvar)
                    if not hasattr(self, 'rotate_mapping'):
                        self.rotate_mapping = {}
                    self.rotate_mapping[xvar] = var

        # Adding variables which may be derived from existing ones
        self.derived_variables = {}
        for m in self.environment_mappings:
            em = self.environment_mappings[m]
            if em['output'][0] not in self.variables and em['input'][0] in self.variables:
                self.logger.debug('Adding variable mapping: %s -> %s' % (em['input'][0], em['output'][0]))
                for v in em['output']:
                    self.variables.append(v)
                    self.derived_variables[v] = em['input']

    def y_is_north(self):
        if self.proj.crs.is_geographic or '+proj=merc' in self.proj.srs:
            return True
        else:
            return False

    def set_buffer_size(self, max_speed, max_vertical_speed=None):
        '''Adjust buffer to minimise data block size needed to cover elements'''
        self.buffer = 0
        pixelsize = self.pixel_size()
        if pixelsize is not None:
            if self.time_step is not None:
                time_step_seconds = self.time_step.total_seconds()
            else:
                time_step_seconds = 3600  # 1 hour if not given
            self.buffer = np.int(np.ceil(max_speed *
                                         time_step_seconds /
                                         pixelsize)) + 2
            self.logger.debug('Setting buffer size %i for reader %s, assuming '
                          'a maximum average speed of %g m/s.' %
                          (self.buffer, self.name, max_speed))

    def pixel_size(self):
        # Find typical pixel size (e.g. for calculating size of buffer)
        if self.projected is True:
            if hasattr(self, 'delta_x'):
                pixelsize = self.delta_x
                if self.proj.crs.is_geographic is True or \
                        ('ob_tran' in self.proj4) or \
                        ('longlat' in self.proj4) or \
                        ('latlon' in self.proj4):
                    pixelsize = pixelsize*111000  # deg to meters
            else:
                pixelsize = None  # Pixel size not defined
        else:
            lons, lats = self.xy2lonlat([self.xmin, self.xmax],
                                        [self.ymin, self.ymin])
            typicalsize = lons
            geod = pyproj.Geod(ellps='WGS84')  # Define an ellipsoid
            dist = geod.inv(lons[0], lats[0],
                            lons[1], lats[1], radians=False)[2]
            pixelsize = dist/self.shape[0]
        return pixelsize

    # TODO:  Duplication from Basemodel, should be unified
    def timer_start(self, category):
        if not hasattr(self, 'timers'):
            self.timers = OrderedDict()
        if not hasattr(self, 'timing'):
            self.timing = OrderedDict()
        if category not in self.timing:
            self.timing[category] = timedelta(0)
        self.timers[category] = datetime.now()

    def timer_end(self, category):
        if self.timers[category] is not None:
            self.timing[category] += datetime.now() - self.timers[category]
        self.timers[category] = None

    def prepare(self, extent, start_time, end_time):
        """Prepare reader for given simulation coverage in time and space."""
        logging.debug('Nothing to prepare for ' + self.name)
        pass  # to be overriden by specific readers

    @abstractmethod
    def get_variables(self, variables, time=None,
                      x=None, y=None, z=None, block=False):
        """Method which must be invoked by any reader (subclass).

        Obtain and return values of the requested variables at all positions
        (x, y, z) closest to given time.

        Arguments:
            variables: string, or list of strings (standard_name) of
                requested variables. These must be provided by reader.
            time: datetime or None, time at which data are requested.
                Can be None (default) if reader/variable has no time
                dimension (e.g. climatology or landmask).
            x, y: float or ndarrays; coordinates of requested points in the
                Spatial Reference System (SRS) of the reader (NB!!)
            z: float or ndarray; vertical position (in meters, positive up)
                of requested points.
                default: 0 m (unless otherwise documented by reader)
            block: bool, see return below

          Returns:
            data: Dictionary
                keywords: variables (string)
                values:
                    - 1D ndarray of len(x) if block=False. Nearest values
                        (neichbour) of requested position are returned.
                    - 3D ndarray encompassing all requested points in
                        x,y,z domain if block=True. It is task of invoking
                        application (OpenDriftSimulation) to perform
                        interpolation in space and time.
        """

    def get_variables_derived(self, variables, *args, **kwargs):
        """Wrapper around get_variables, adding derived"""
        if isinstance(variables, basestring):
            variables = [variables]
        if not isinstance(variables, list):
            variables = list(variables)
        derive_variables = False
        for var in variables:
            if var in self.derived_variables:
                fromvars = self.derived_variables[var]
                for v in fromvars:
                    variables.append(v)
                # Removing the derived variable name
                variables = [v for v in variables if v != var]
                derive_variables = True

        env = self.get_variables(variables, *args, **kwargs)
        if derive_variables is True:
            self.calculate_derived_environment_variables(env)

        return env

    def _get_variables(self, variables, profiles, profiles_depth,
                       time, x, y, z, block):
        """Wrapper around reader-specific function get_variables()

        Performs some common operations which should not be duplicated:
        - monitor time spent by this reader
        - convert any numpy arrays to masked arrays
        """
        self.logger.debug('Fetching variables from ' + self.name)
        self.timer_start('reading')

        if profiles is not None and block is True:
            # If profiles are requested for any parameters, we
            # add two fake points at the end of array to make sure that the
            # requested block has the depth range required for profiles
            x = np.append(x, [x[-1], x[-1]])
            y = np.append(y, [y[-1], y[-1]])
            z = np.append(z, [profiles_depth[0], profiles_depth[1]])
        env = self.get_variables_derived(variables, time, x, y, z, block)

        # Make sure x and y are floats (and not e.g. int64)
        if 'x' in env.keys():
            env['x'] = np.array(env['x'], dtype=np.float)
            env['y'] = np.array(env['y'], dtype=np.float)

        for variable in variables:
            # Convert any masked arrays to NumPy arrays
            if isinstance(env[variable], np.ma.MaskedArray):
                env[variable] = env[variable].filled(np.nan)
            # Mask values outside valid_min, valid_max (self.standard_names)
            if variable in standard_names.keys():
                if isinstance(env[variable], list):
                    self.logger.warning('Skipping min-max checking for ensemble data')
                    continue
                with np.errstate(invalid='ignore'):
                    invalid_indices = np.logical_and(
                        np.isfinite(env[variable]), np.logical_or(
                        env[variable]<standard_names[variable]['valid_min'],
                        env[variable]>standard_names[variable]['valid_max']))
                if np.sum(invalid_indices) > 0:
                    invalid_values = env[variable][invalid_indices]
                    self.logger.warning('Invalid values (%s to %s) found for %s, replacing with NaN' % (invalid_values.min(), invalid_values.max(), variable))
                    self.logger.warning('(allowed range: [%s, %s])' %
                                    (standard_names[variable]['valid_min'],
                                     standard_names[variable]['valid_max']))
                    env[variable][invalid_indices] = np.nan

        # Convolve arrays with a kernel, if reader.convolve is set
        if hasattr(self, 'convolve'):
            from scipy import ndimage
            N = self.convolve
            if isinstance(N, (int, np.integer)):
                kernel = np.ones((N, N))
                kernel = kernel/kernel.sum()
            else:
                kernel = N
            self.logger.debug('Convolving variables with kernel: %s' % kernel)
            for variable in env.keys():
                if variable in ['x', 'y', 'z', 'time']:
                    pass
                else:
                    if env[variable].ndim == 2:
                        env[variable] = ndimage.convolve(
                            env[variable], kernel, mode='nearest')
                    elif env[variable].ndim == 3:
                        env[variable] = ndimage.convolve(
                            env[variable], kernel[:,:,None],
                            mode='nearest')

        self.timer_end('reading')

        return env

    environment_mappers = []

    def wind_from_speed_and_direction(self, env):
        north_wind = env['wind_speed']*np.cos(
            np.radians(env['wind_to_direction']))
        east_wind = env['wind_speed']*np.sin(
            np.radians(env['wind_to_direction']))
        env['x_wind'] = east_wind
        env['y_wind'] = north_wind
        # Rotating might be necessary generally
        #x,y = np.meshgrid(env['x'], env['y'])
        #env['x_wind'], env['y_wind'] = self.rotate_vectors(
        #    x, y,
        #    east_wind, north_wind,
        #    None, self.proj)

    environment_mappings = {
        'wind_from_speed_and_direction': {
            'input': ['wind_speed', 'wind_to_direction'],
            'output': ['x_wind', 'y_wind'],
            'method': wind_from_speed_and_direction},
        'testvar': {
            'input': ['sea_ice_thickness'],
            'output': ['istjukkleik']}
        }

    def calculate_derived_environment_variables(self, env):

        if 'x_wind' in self.derived_variables and 'wind_speed' in env.keys():
            self.wind_from_speed_and_direction(env)

    def get_variables_interpolated(self, variables, profiles=None,
                                   profiles_depth=None, time=None,
                                   lon=None, lat=None, z=None,
                                   block=False, rotate_to_proj=None):

        self.timer_start('total')
        self.timer_start('preparing')
        # Raise error if time not not within coverage of reader
        if not self.covers_time(time):
            raise ValueError('%s is outside time coverage (%s - %s) of %s' %
                             (time, self.start_time, self.end_time, self.name))

        # Check which particles are covered (indep of time)
        ind_covered, reader_x, reader_y = self.covers_positions(lon, lat, z)
        if len(ind_covered) == 0:
            raise ValueError(('All %s particles (%.2f-%.2fE, %.2f-%.2fN) ' +
                              'are outside domain of %s (%s)') %
                             (len(lon), lon.min(), lon.max(), lat.min(),
                              lat.max(), self.name, self.coverage_string()))

        # Find reader time_before/time_after
        time_nearest, time_before, time_after, i1, i2, i3 = \
            self.nearest_time(time)
        self.logger.debug('Reader time:\n\t\t%s (before)\n\t\t%s (after)' %
                      (time_before, time_after))
        # For variables which are not time dependent, we do not care about time
        static_variables = ['sea_floor_depth_below_sea_level', 'land_binary_mask']
        if time == time_before or all(v in static_variables for v in variables):
            time_after = None

        if len(z) == 1 and len(lon) > 1:
            z = z.copy()*np.ones(lon.shape)
        z = z.copy()[ind_covered]  # Send values and not reference
                                   # to avoid modifications
        if block is False or self.return_block is False:
            # Analytical reader, continous in space and time
            self.timer_end('preparing')
            env_before = self._get_variables(variables, profiles,
                                             profiles_depth,
                                             time,
                                             #time_before,
                                             reader_x, reader_y, z,
                                             block=block)
            self.logger.debug('Fetched env-before')
            self.timer_start('preparing')

        else:
            block_before = block_after = None
            blockvariables_before = variables
            blockvars_before = str(variables)
            blockvariables_after = variables
            blockvars_after = str(variables)
            for blockvars in self.var_block_before:
                if all(v in blockvars for v in variables):
                    block_before = self.var_block_before[blockvars]
                    blockvariables_before = block_before.data_dict.keys()
                    blockvars_before = blockvars
                    break
                blockvariables_before = variables
                blockvars_before = str(variables)
            for blockvars in self.var_block_after:
                if all(v in blockvars for v in variables):
                    block_after = self.var_block_after[blockvars]
                    blockvariables_after = block_after.data_dict.keys()
                    blockvars_after = blockvars
                    break

            # Swap before- and after-blocks if matching times
            if block_before is not None and block_after is not None:
                if block_before.time != time_before:
                    if block_after.time == time_before:
                        block_before = block_after
                        self.var_block_before[blockvars_before] = block_before
                if block_after.time != time_after:
                    if block_before.time == time_before:
                        block_after = block_before
                        self.var_block_after[blockvars_after] = block_after
            # Fetch data, if no buffer is available
            if block_before is None or \
                    block_before.time != time_before:
                self.timer_end('preparing')
                reader_data_dict = \
                    self._get_variables(blockvariables_before, profiles,
                                        profiles_depth, time_before,
                                        reader_x, reader_y, z,
                                        block=block)
                self.timer_start('preparing')
                self.var_block_before[blockvars_before] = \
                    ReaderBlock(reader_data_dict,
                                interpolation_horizontal=self.interpolation)
                try:
                    len_z = len(self.var_block_before[blockvars_before].z)
                except:
                    len_z = 1
                self.logger.debug(('Fetched env-block (size %ix%ix%i) ' +
                              'for time before (%s)') %
                              (len(self.var_block_before[blockvars_before].x),
                               len(self.var_block_before[blockvars_before].y),
                               len_z, time_before))
                block_before = self.var_block_before[blockvars_before]
            if block_after is None or block_after.time != time_after:
                if time_after is None:
                    self.var_block_after[blockvars_after] = \
                        block_before
                else:
                    self.timer_end('preparing')
                    reader_data_dict = \
                        self._get_variables(blockvariables_after, profiles,
                                            profiles_depth, time_after,
                                            reader_x, reader_y, z,
                                            block=block)
                    self.timer_start('preparing')
                    self.var_block_after[blockvars_after] = \
                        ReaderBlock(
                            reader_data_dict,
                            interpolation_horizontal=self.interpolation)
                    try:
                        len_z = len(self.var_block_after[blockvars_after].z)
                    except:
                        len_z = 1

                    self.logger.debug(('Fetched env-block (size %ix%ix%i) ' +
                                  'for time after (%s)') %
                                  (len(self.var_block_after[blockvars_after].x),
                                   len(self.var_block_after[blockvars_after].y),
                                   len_z, time_after))
                    block_after = self.var_block_after[blockvars_after]

            if (block_before is not None and block_before.covers_positions(
                reader_x, reader_y) is False) or (\
                block_after is not None and block_after.covers_positions(
                    reader_x, reader_y) is False):
                self.logger.warning('Data block from %s not large enough to '
                                'cover element positions within timestep. '
                                'Buffer size (%s) must be increased.' %
                                (self.name, str(self.buffer)))

            self.timer_end('preparing')
            ############################################################
            # Interpolate before/after blocks onto particles in space
            ############################################################
            self.timer_start('interpolation')
            self.logger.debug('Interpolating before (%s) in space  (%s)' %
                          (block_before.time, self.interpolation))
            env_before, env_profiles_before = block_before.interpolate(
                    reader_x, reader_y, z, variables,
                    profiles, profiles_depth)

            if (time_after is not None) and (time_before != time):
                self.logger.debug('Interpolating after (%s) in space  (%s)' %
                              (block_after.time, self.interpolation))
                env_after, env_profiles_after = block_after.interpolate(
                        reader_x, reader_y, z, variables,
                        profiles, profiles_depth)

            self.timer_end('interpolation')

        #######################
        # Time interpolation
        #######################
        self.timer_start('interpolation_time')
        env_profiles = None
        if (time_after is not None) and (time_before != time) and self.return_block is True:
            weight_after = ((time - time_before).total_seconds() /
                            (time_after - time_before).total_seconds())
            self.logger.debug(('Interpolating before (%s, weight %.2f) and'
                           '\n\t\t      after (%s, weight %.2f) in time') %
                          (block_before.time, 1 - weight_after,
                           block_after.time, weight_after))
            env = {}
            for var in variables:
                # Weighting together, and masking invalid entries
                env[var] = np.ma.masked_invalid((env_before[var] *
                                                (1 - weight_after) +
                                                env_after[var] * weight_after))
            # Interpolating vertical profiles in time
            if profiles is not None:
                env_profiles = {}
                self.logger.debug('Interpolating profiles in time')
                # Truncating layers not present both before and after
                numlayers = np.minimum(len(env_profiles_before['z']),
                                       len(env_profiles_after['z']))
                env_profiles['z'] = env_profiles_before['z'][0:numlayers+1]
                for var in env_profiles_before.keys():
                    if var == 'z':
                        continue
                    env_profiles_before[var]=np.atleast_2d(env_profiles_before[var])
                    env_profiles_after[var]=np.atleast_2d(env_profiles_after[var])
                    env_profiles[var] = (
                        env_profiles_before[var][0:numlayers, :] *
                        (1 - weight_after) +
                        env_profiles_after[var][0:numlayers, :]*weight_after)
            else:
                env_profiles = None

        else:
            self.logger.debug('No time interpolation needed - right on time.')
            env = env_before
            if profiles is not None:
                if 'env_profiles_before' in locals():
                    env_profiles = env_profiles_before
                else:
                    # Copying data from environment to vertical profiles
                    env_profiles = {'z': profiles_depth}
                    for var in profiles:
                        env_profiles[var] = np.ma.array([env[var], env[var]])
        self.timer_end('interpolation_time')

        ####################
        # Rotate vectors
        ####################
        if rotate_to_proj is not None:
            if self.simulation_SRS is True:
                self.logger.debug('Reader SRS is the same as calculation SRS - '
                              'rotation of vectors is not needed.')
            else:
                vector_pairs = []
                for var in variables:
                    for vector_pair in vector_pairs_xy:
                        if var in vector_pair:
                            counterpart = list(set(vector_pair) -
                                               set([var]))[0]
                            if counterpart in variables:
                                vector_pairs.append(vector_pair)
                            else:
                                self.logger.warning('Missing component of vector pair, cannot rotate:' +
                                         counterpart)
                # Extract unique vector pairs
                vector_pairs = [list(x) for x in set(tuple(x)
                                for x in vector_pairs)]

                if len(vector_pairs) > 0:
                    self.timer_start('rotating vectors')
                    for vector_pair in vector_pairs:
                        env[vector_pair[0]], env[vector_pair[1]] = \
                            self.rotate_vectors(reader_x, reader_y,
                                                env[vector_pair[0]],
                                                env[vector_pair[1]],
                                                self.proj, rotate_to_proj)
                        if profiles is not None and vector_pair[0] in profiles:
                            env_profiles[vector_pair[0]], env_profiles[vector_pair[1]] = \
                                    self.rotate_vectors (reader_x, reader_y,
                                            env_profiles[vector_pair[0]],
                                            env_profiles[vector_pair[1]],
                                            self.proj,
                                            rotate_to_proj)


                    self.timer_end('rotating vectors')

        # Masking non-covered pixels
        self.timer_start('masking')
        if len(ind_covered) != len(lon):
            self.logger.debug('Masking %i elements outside coverage' %
                          (len(lon)-len(ind_covered)))
            for var in variables:
                tmp = np.nan*np.ones(lon.shape)
                tmp[ind_covered] = env[var].copy()
                env[var] = np.ma.masked_invalid(tmp)
                # Filling also fin missing columns
                # for env_profiles outside coverage
                if env_profiles is not None and var in env_profiles.keys():
                    tmp = np.nan*np.ones((env_profiles[var].shape[0],
                                          len(lon)))
                    tmp[:, ind_covered] = env_profiles[var].copy()
                    env_profiles[var] = np.ma.masked_invalid(tmp)

        self.timer_end('masking')
        self.timer_end('total')
        return env, env_profiles

    def rotate_variable_dict(self, variables, proj_from='+proj=latlong', proj_to=None):
        for vectorpair in vector_pairs_xy:
            if vectorpair[0] in self.rotate_mapping and vectorpair[0] in variables.keys():
                if proj_to is None:
                    proj_to = self.proj
                self.logger.debug('Rotating vector from east/north to xy orientation: ' + str(vectorpair))
                variables[vectorpair[0]], variables[vectorpair[1]] = self.rotate_vectors(
                    variables['x'], variables['y'],
                    variables[vectorpair[0]], variables[vectorpair[1]],
                    proj_from, proj_to)

    def rotate_vectors(self, reader_x, reader_y,
                       u_component, v_component,
                       proj_from, proj_to):
        """Rotate vectors from one srs to another."""

        if type(proj_from) is str:
            proj_from = pyproj.Proj(proj_from)
        if type(proj_from) is not pyproj.Proj:
            proj_from = pyproj.Proj('+proj=latlong +R=6370997.0 +ellps=WGS84')
            reader_x, reader_y = self.xy2lonlat(reader_x, reader_y)
        if type(proj_to) is str:
            proj_to = pyproj.Proj(proj_to)

        if proj_from.crs.is_geographic:
                delta_y = .1  # 0.1 degree northwards
        else:
            delta_y = 10  # 10 m along y-axis
        
        transformer = pyproj.Transformer.from_proj(proj_from, proj_to)
        x2, y2 = transformer.transform(reader_x, reader_y)
        x2_delta, y2_delta = transformer.transform(reader_x, reader_y + delta_y)

        if proj_to.crs.is_geographic:
            geod = pyproj.Geod(ellps='WGS84')
            rot_angle_vectors_rad = np.radians(
                geod.inv(x2, y2, x2_delta, y2_delta)[0])
        else:
            rot_angle_vectors_rad = np.arctan2(x2_delta - x2, y2_delta - y2)
        self.logger.debug('Rotating vectors between %s and %s degrees.' %
                      (np.degrees(rot_angle_vectors_rad).min(),
                       np.degrees(rot_angle_vectors_rad).max()))
        rot_angle_rad = - rot_angle_vectors_rad
        u_rot = (u_component*np.cos(rot_angle_rad) -
                 v_component*np.sin(rot_angle_rad))
        v_rot = (u_component*np.sin(rot_angle_rad) +
                 v_component*np.cos(rot_angle_rad))

        return u_rot, v_rot

    def xy2lonlat(self, x, y):
        """Calculate x,y in own projection from given lon,lat (scalars/arrays).
        """
        if self.projected is True:
            if self.proj.crs.is_geographic:
                if 'ob_tran' in self.proj4:
                    self.logger.debug('NB: Converting degrees to radians ' +
                                 'due to ob_tran srs')
                    x = np.radians(np.array(x))
                    y = np.radians(np.array(y))
                    return self.proj(x, y, inverse=True)
                else:
                    return x, y
            else:
                return self.proj(x, y, inverse=True)
        else:
            np.seterr(invalid='ignore')  # Disable warnings for nan-values
            y = np.atleast_1d(np.array(y))
            x = np.atleast_1d(np.array(x))

            # NB: mask coordinates outside domain
            x[x < self.xmin] = np.nan
            x[x > self.xmax] = np.nan
            y[y < self.ymin] = np.nan
            y[y < self.ymin] = np.nan

            lon = map_coordinates(self.lon, [y, x], order=1,
                                  cval=np.nan, mode='nearest')
            lat = map_coordinates(self.lat, [y, x], order=1,
                                  cval=np.nan, mode='nearest')
            return (lon, lat)

    def lonlat2xy(self, lon, lat):
        """Calculate lon,lat from given x,y (scalars/arrays) in own projection.
        """
        
        # Two methods needed for parallelisation
        def get_x(lon_part, lat_part, num):
            out_x[num] = self.spl_x(lon_part, lat_part)
        def get_y(lon_part, lat_part, num):
            out_y[num] = self.spl_y(lon_part, lat_part)

        if self.projected is True:
            if 'ob_tran' in self.proj4:
                x, y = self.proj(lon, lat, inverse=False)
                return np.degrees(x), np.degrees(y)
            elif self.proj.crs.is_geographic:
                return lon, lat
            else:
                x, y = self.proj(lon, lat, inverse=False)
                return x, y
        else:
            # For larger arrays, we split and calculate in parallel
            # The number of CPUs to use can be improved/optimised
            num_elements = len(np.atleast_1d(lon))
            if num_elements > 100 and not hasattr(self, 'multiprocessing_fail'):
                try:
                    nproc = 2
                    if num_elements > 1000:
                        nproc = 16
                    if num_elements > 100000:
                        nproc = 32
                    if num_elements > 1000000:
                        nproc = 64
                    cpus = cpu_count()
                    nproc = np.minimum(nproc, cpus-1)
                    nproc = np.maximum(2, nproc)
                    self.logger.debug('Running lonlat2xy in parallel, using %i of %i CPUs'
                                  % (nproc, cpus))
                    split_lon = np.array_split(lon, nproc)
                    split_lat = np.array_split(lat, nproc)
                    out_x = Manager().dict()
                    out_y = Manager().dict()
                    processes = []
                    for i in range(nproc):
                        processes.append(Process(target=get_x, args=
                                         (split_lon[i], split_lat[i], i)))
                        processes.append(Process(target=get_y, args=
                                         (split_lon[i], split_lat[i], i)))
                    [p.start() for p in processes]
                    [p.join() for p in processes]
                    x = np.concatenate(out_x)
                    y = np.concatenate(out_y)
                    self.logger.debug('Completed lonlat2xy in parallel')
                    return (x, y)
                except Exception as e:
                    self.logger.warning('Parallelprocessing failed:')
                    self.logger.warning(e)
                    self.multiprocessing_fail = True
            else:
                if hasattr(self, 'multiprocessing_fail'):
                    self.logger.warning('Multiprocessing has previously failed, reverting to using single processor for lonlat -> xy')
                else:
                    # For smaller arrays, we run sequentially
                    pass
            self.logger.debug('Calculating lonlat->xy sequentially')
            x = self.spl_x(lon, lat)
            y = self.spl_y(lon, lat)
            return (x, y)

    def y_azimuth(self, lon, lat):
        """Calculate azimuth orientation of the y-axis of the reader SRS."""
        x0, y0 = self.lonlat2xy(lon, lat)
        distance = 1000.0  # Create points 1 km away to determine azimuth
        lon_2, lat_2 = self.xy2lonlat(x0, y0 + distance)
        geod = pyproj.Geod(ellps='WGS84')  # Define an ellipsoid
        y_az = geod.inv(lon_2, lat_2, lon, lat, radians=False)
        dist = y_az[2]
        return y_az[0]

    def covers_time(self, time):

        if self.always_valid is True:
            return True
        if self.start_time is None:
            return True  # No time limitations of reader
        if (time < self.start_time) or (time > self.end_time):
            return False
        else:
            return True

    def covers_positions(self, lon, lat, z=0):
        """Return indices of input points covered by reader."""

        z = np.atleast_1d(z)
        # Calculate x,y coordinates from lon,lat
        x, y = self.lonlat2xy(lon, lat)

        # Only checking vertical coverage if zmin, zmax is defined
        zmin = -np.inf
        zmax = np.inf
        if hasattr(self, 'zmin') and self.zmin is not None:
            zmin = self.zmin
        if hasattr(self, 'zmax') and self.zmax is not None:
            zmax = self.zmax

        if self.global_coverage():
            # We need only check north-south and z coverage
            indices = np.where((y >= self.ymin) & (y <= self.ymax) &
                               (z >= zmin) & (z <= zmax))[0]
        else:
            indices = np.where((x >= self.xmin) & (x <= self.xmax) &
                               (y >= self.ymin) & (y <= self.ymax) &
                               (z >= zmin) & (z <= zmax))[0]

        try:
            return indices, x[indices], y[indices]
        except:
            return indices, x, y

    def global_coverage(self):
        """Return True if global coverage east-west"""

        if self.proj.crs.is_geographic is True and hasattr(self, 'delta_x'):
            if (self.xmin - self.delta_x <= 0) and (
                self.xmax + self.delta_x >= 360):
                return True  # Global 0 to 360
            if (self.xmin - self.delta_x <= -180) and (
                self.xmax + self.delta_x >= 180):
                return True  # Global -180 to 180

        return False

    def coverage_string(self):
        """Coverage of reader to be reported as string for debug output"""
        corners = self.xy2lonlat([self.xmin, self.xmin, self.xmax, self.xmax],
                                 [self.ymax, self.ymin, self.ymax, self.ymin])
        return '%.2f-%.2fE, %.2f-%.2fN' % (
                    np.min(corners[0]), np.max(corners[0]),
                    np.min(corners[1]), np.max(corners[1]))

    def check_arguments(self, variables, time, x, y, z):
        """Check validity of arguments input to method get_variables.

        Checks that requested positions and time are within coverage of
        this reader, and that it can provide the requested variable(s).
        Returns the input arguments, possibly modified/corrected (below)

        Arguments:
            See function get_variables for definition.

        Returns:
            variables: same as input, but converted to list if given as string.
            time: same as input, or start_time of reader if given as None.
            x, y, z: same as input, but converted to ndarrays
                if given as scalars.
            outside: boolean array which is True for any particles outside
                the spatial domain covered by this reader.

        Raises:
            ValueError:
                - if requested time is outside coverage of reader.
                - if all requested positions are outside coverage of reader.
        """

        # Check time
        if time is None:
            time = self.start_time  # Use first timestep, if not given

        # Convert variables to list and x,y to ndarrays
        if isinstance(variables, basestring):
            variables = [variables]
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        if z is not None:
            z = np.asarray(z)

        for variable in variables:
            if variable not in self.variables:
                raise ValueError('Variable not available: ' + variable +
                                 '\nAvailable parameters are: ' +
                                 str(self.variables))
        if (self.start_time is not None and time < self.start_time) and self.always_valid is False:
            raise ValueError('Requested time (%s) is before first available '
                             'time (%s) of %s' % (time, self.start_time,
                                                  self.name))
        if (self.end_time is not None and time > self.end_time) and self.always_valid is False:
            raise ValueError('Requested time (%s) is after last available '
                             'time (%s) of %s' % (time, self.end_time,
                                                  self.name))
        if self.global_coverage():
            outside = np.where(~np.isfinite(x+y) |
                               (y < self.ymin) | (y > self.ymax))[0]
        else:
            outside = np.where(~np.isfinite(x+y) |
                               (x < self.xmin) | (x > self.xmax) |
                               (y < self.ymin) | (y > self.ymax))[0]
        if np.size(outside) == np.size(x):
            lon, lat = self.xy2lonlat(x, y)
            raise ValueError(('Argcheck: all %s particles (%.2f-%.2fE, ' +
                              '%.2f-%.2fN) are outside domain of %s (%s)') %
                             (len(lon), lon.min(), lon.max(), lat.min(),
                              lat.max(), self.name, self.coverage_string()))

        return variables, time, x, y, z, outside

    def nearest_time(self, time):
        """Return nearest times before and after the requested time.

        Returns:
            nearest_time: datetime
            time_before: datetime
            time_after: datetime
            indx_nearest: int
            indx_before: int
            indx_after: int
        """
        if self.start_time == self.end_time:
            return self.start_time, self.start_time, None, 0, 0, 0
        if self.start_time is None:
            return None, None, None, None, None, None
        if hasattr(self, 'times'):  # Time as array, possibly with holes
            indx_before = np.max((0, bisect_left(self.times, time) - 1))
            if self.times[indx_before + 1] == time:
                # Correction needed when requested time exists in times
                indx_before = indx_before + 1
            time_before = self.times[indx_before]
            indx_after = np.minimum(indx_before + 1,
                                    len(self.times) - 1)  # At the end
            time_after = self.times[indx_after]
            if (time - time_before) < (time_after - time):
                indx_nearest = indx_before
            else:
                indx_nearest = indx_after
            nearest_time = self.times[indx_nearest]
        else:  # Time step is constant (no holes)
            if self.time_step is None:
                return None, None, None, None, None, None
            indx = float((time - self.start_time).total_seconds()) / \
                float(self.time_step.total_seconds())
            indx_nearest = int(round(indx))
            nearest_time = self.start_time + indx_nearest*self.time_step
            indx_before = int(np.floor(indx))
            time_before = self.start_time + indx_before*self.time_step
            indx_after = int(np.ceil(indx))
            time_after = self.start_time + indx_after*self.time_step
        return nearest_time, time_before, time_after,\
            indx_nearest, indx_before, indx_after

    def index_of_closest_z(self, requested_z):
        """Return (internal) index of z closest to requested z.

        Thickness of layers (of ocean model) are not assumed to be constant.
        """
        ind_z = [np.abs(np.subtract.outer(
            self.z, requested_z)).argmin(0)]
        return ind_z, self.z[ind_z]

    def indices_min_max_z(self, z):
        """
        Return min and max indices of internal vertical dimension,
        covering the requested vertical positions.
        Needed when block is requested (True).

        Arguments:
            z: ndarray of floats, in meters
        """
        minIndex = (self.z <= z.min()).argmin() - 1
        maxIndex = (self.z >= z.max()).argmax()
        return minIndex, maxIndex

    def domain_grid(self, npoints=1000):
        """Return arrays of lon,lat points spread evenly over reader domain."""
        numx = np.floor(np.sqrt(npoints))
        numy = np.floor(np.sqrt(npoints))
        x = np.linspace(self.xmin, self.xmax - self.delta_x, numx)
        y = np.linspace(self.ymin, self.ymax - self.delta_y, numy)
        x, y = np.meshgrid(x, y)
        lons, lats = self.xy2lonlat(x, y)
        return lons, lats

    def __repr__(self):
        """String representation of the current reader."""
        outStr = '===========================\n'
        outStr += 'Reader: ' + self.name + '\n'
        outStr += 'Projection: \n  ' + self.proj4 + '\n'
        if self.proj.crs.is_geographic:
            if self.projected is False:
                outStr += 'Coverage: [pixels]\n'
            else:
                outStr += 'Coverage: [degrees]\n'
        else:
            if self.projected is False:
                outStr += 'Coverage: [pixels]\n'
            else:
                outStr += 'Coverage: [m]\n'
        shape = self.shape
        if shape is None:
            outStr += '  xmin: %f   xmax: %f\n' % (self.xmin, self.xmax)
            outStr += '  ymin: %f   ymax: %f\n' % (self.ymin, self.ymax)
        else:
            outStr += '  xmin: %f   xmax: %f   step: %g   numx: %i\n' % \
                (self.xmin, self.xmax, self.delta_x or 0, shape[0])
            outStr += '  ymin: %f   ymax: %f   step: %g   numy: %i\n' % \
                (self.ymin, self.ymax, self.delta_y or 0, shape[1])
        corners = self.xy2lonlat([self.xmin, self.xmin, self.xmax, self.xmax],
                                 [self.ymax, self.ymin, self.ymax, self.ymin])
        outStr += '  Corners (lon, lat):\n'
        outStr += '    (%6.2f, %6.2f)  (%6.2f, %6.2f)\n' % \
            (corners[0][0],
             corners[1][0],
             corners[0][2],
             corners[1][2])
        outStr += '    (%6.2f, %6.2f)  (%6.2f, %6.2f)\n' % \
            (corners[0][1],
             corners[1][1],
             corners[0][3],
             corners[1][3])
        if hasattr(self, 'z'):
            np.set_printoptions(suppress=True)
            outStr += 'Vertical levels [m]: \n  ' + str(self.z) + '\n'
        elif hasattr(self, 'sigma'):
            outStr += 'Vertical levels [sigma]: \n  ' + str(self.sigma) + '\n'
        else:
            outStr += 'Vertical levels [m]: \n  Not specified\n'
        outStr += 'Available time range:\n'
        outStr += '  start: ' + str(self.start_time) + \
                  '   end: ' + str(self.end_time) + \
                  '   step: ' + str(self.time_step) + '\n'
        if self.start_time is not None and self.time_step is not None:
            outStr += '    %i times (%i missing)\n' % (
                      self.expected_time_steps, self.missing_time_steps)
        if hasattr(self, 'realizations'):
            outStr += 'Variables (%i ensemble members):\n' % len(self.realizations)
        else:
            outStr += 'Variables:\n'
        for variable in self.variables:
            if variable in self.derived_variables:
                outStr += '  ' + variable + ' - derived from ' + \
                    str(self.derived_variables[variable]) + '\n'
            else:
                outStr += '  ' + variable + '\n'
        outStr += '===========================\n'
        outStr += self.performance()

        return outStr

    def performance(self):
        '''Report the time spent on various tasks'''
        outStr = ''
        if hasattr(self, 'timing'):
            for cat, time in self.timing.items():
                time = str(time)[0:str(time).find('.') + 2]
                outStr += '%10s  %s\n' % (time, cat)
        return outStr

    def clip_boundary_pixels(self, numpix):
        '''Trim some (potentially bad) pixels along boundary'''
        self.logger.info('Trimming %i pixels from boundary' % numpix)
        self.xmin = self.xmin+numpix*self.delta_x
        self.xmax = self.xmax-numpix*self.delta_x
        self.ymin = self.ymin+numpix*self.delta_y
        self.ymax = self.ymax-numpix*self.delta_y
        self.shape = tuple([self.shape[0]-2*numpix,
                            self.shape[1]-2*numpix])
        self.clipped = numpix

    def plot(self, variable=None, vmin=None, vmax=None,
             filename=None, title=None, buffer=1, lscale='auto'):
        """Plot geographical coverage of reader."""

        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        from opendrift_landmask_data import Landmask
        fig = plt.figure()

        corners = self.xy2lonlat([self.xmin, self.xmin, self.xmax, self.xmax],
                                 [self.ymax, self.ymin, self.ymax, self.ymin])
        lonmin = np.min(corners[0]) - buffer*2
        lonmax = np.max(corners[0]) + buffer*2
        latmin = np.min(corners[1]) - buffer
        latmax = np.max(corners[1]) + buffer
        latspan = latmax - latmin

        # Initialise map
        if latspan < 90:
            # Stereographic projection centred on domain, if small domain
            x0 = (self.xmin + self.xmax) / 2
            y0 = (self.ymin + self.ymax) / 2
            lon0, lat0 = self.xy2lonlat(x0, y0)
            sp = ccrs.Stereographic(central_longitude=lon0, central_latitude=lat0)
            ax = fig.add_subplot(1, 1, 1, projection=sp)
            corners_stere = sp.transform_points(ccrs.PlateCarree(), np.array(corners[0]), np.array(corners[1]))
        else:
            # Global map if reader domain is large
            sp = ccrs.Mercator()
            ax = fig.add_subplot(1, 1, 1, projection=sp)

        # GSHHS coastlines
        f = cfeature.GSHHSFeature(scale=lscale, levels=[1],
                                  facecolor=cfeature.COLORS['land'])
        ax.add_geometries(
            f.intersecting_geometries([lonmin, lonmax, latmin, latmax]),
            ccrs.PlateCarree(),
            facecolor=cfeature.COLORS['land'],
            edgecolor='black')

        gl = ax.gridlines(ccrs.PlateCarree())
        gl.top_labels = False

        # Get boundary
        npoints = 10  # points per side
        x = np.array([])
        y = np.array([])
        x = np.concatenate((x, np.linspace(self.xmin, self.xmax, npoints)))
        y = np.concatenate((y, [self.ymin]*npoints))
        x = np.concatenate((x, [self.xmax]*npoints))
        y = np.concatenate((y, np.linspace(self.ymin, self.ymax, npoints)))
        x = np.concatenate((x, np.linspace(self.xmax, self.xmin, npoints)))
        y = np.concatenate((y, [self.ymax]*npoints))
        x = np.concatenate((x, [self.xmin]*npoints))
        y = np.concatenate((y, np.linspace(self.ymax, self.ymin, npoints)))
        # from x/y vectors create a Patch to be added to map
        lon, lat = self.xy2lonlat(x, y)
        lat[lat>89] = 89.
        lat[lat<-89] = -89.
        p = sp.transform_points(ccrs.PlateCarree(), lon, lat)
        xsp = p[:, 0]
        ysp = p[:, 1]

        if variable is None:
            boundary = Polygon(list(zip(xsp, ysp)), alpha=0.5, ec='k', fc='b',
                               zorder=100)
            ax.add_patch(boundary)
            buf = (xsp.max()-xsp.min())*.1  # Some whitespace around polygon
            buf = 0
            try:
                ax.set_extent([xsp.min()-buf, xsp.max()+buf, ysp.min()-buf, ysp.max()+buf], crs=sp)
            except:
                pass
        if title is None:
            plt.title(self.name)
        else:
            plt.title(title)
        plt.xlabel('Time coverage: %s to %s' %
                   (self.start_time, self.end_time))

        if variable is not None:
            rx = np.array([self.xmin, self.xmax])
            ry = np.array([self.ymin, self.ymax])
            data = self.get_variables_derived(variable, self.start_time,
                                      rx, ry, block=True)
            rx, ry = np.meshgrid(data['x'], data['y'])
            rlon, rlat = self.xy2lonlat(rx, ry)
            data[variable] = np.ma.masked_invalid(data[variable])
            if hasattr(self, 'convolve'):
                from scipy import ndimage
                N = self.convolve
                if isinstance(N, (int, np.integer)):
                    kernel = np.ones((N, N))
                    kernel = kernel/kernel.sum()
                else:
                    kernel = N
                self.logger.debug('Convolving variables with kernel: %s' % kernel)
                data[variable] = ndimage.convolve(
                            data[variable], kernel, mode='nearest')
            if data[variable].ndim > 2:
                self.logger.warning('Ensemble data, plotting only first member')
                data[variable] = data[variable][0,:,:]
            mappable = ax.pcolormesh(rlon, rlat, data[variable], vmin=vmin, vmax=vmax,
                                     transform=ccrs.PlateCarree())
            cbar = fig.colorbar(mappable, orientation='horizontal', pad=.05, aspect=30, shrink=.4)
            cbar.set_label(variable)

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except:
            pass

        if filename is not None:
            plt.savefig(filename)
            plt.close()
        else:
            plt.show()

    def get_timeseries_at_position(self, lon, lat, variables=None,
                                   start_time=None, end_time=None, times=None):
        """ Get timeseries of variables from this reader at given position.
        """

        if times is None:
            if start_time is None:
                start_time = self.start_time
            if end_time is None:
                end_time = self.end_time
            times = [t for t in self.times if t >= start_time and t<= end_time]

        if variables is None:
            variables = self.variables

        if len(self.covers_positions(lon=lon, lat=lat)[0]) == 0:
            return None

        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)
        if len(lon) == 1:
            lon = lon[0]*np.ones(len(times))
            lat = lat[0]*np.ones(len(times))

        data = {'time': times}
        for var in variables:
            data[var] = np.zeros(len(times))

        for i, time in enumerate(times):
            closest_time = min(self.times, key=lambda d: abs(d - time))
            print(time, closest_time, 'Time, Closest time')
            d = self.get_variables_interpolated(
                lon=np.atleast_1d(lon[i]), lat=np.atleast_1d(lat[i]), z=np.atleast_1d(0),
                time=closest_time, variables=variables, rotate_to_proj='+proj=latlong')[0]
            for var in variables:
                data[var][i] = d[var][0]

        return(data)
