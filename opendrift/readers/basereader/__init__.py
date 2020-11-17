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
from multiprocessing import Process, Manager, cpu_count

from scipy.interpolate import LinearNDInterpolator
import numpy as np
import pyproj

from .structured import StructuredReader
from .unstructured import UnstructuredReader
from .continuous import ContinuousReader
from .variables import Variables
from .fakeproj import fakeproj
from .consts import *

from opendrift.readers.interpolation import ReaderBlock

class BaseReader(Variables):
    """
    An abstract reader. Implementors provide a method to read data and specify how it is interpolated.
    """

    __metaclass__ = ABCMeta

    return_block = True  # By default, all readers should be
                         # capable of returning blocks of data

    # Default interpolation method, see function interpolate_block()
    interpolation = 'linearNDFast'

    convolve = None # Convolution kernel or kernel size

    verticalbuffer = 1  # To be overridden by application as needed

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
        'salinity_vertical_diffusion_coefficient' : 'ocean_vertical_diffusivity',
        'sea_floor_depth_below_geoid' : 'sea_floor_depth_below_sea_level'
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
        """Common constructor for all readers"""
        super().__init__()

        # Dictionaries to store blocks of data for reuse (buffering)
        self.var_block_before = {}  # Data for last timestep before present
        self.var_block_after = {}   # Data for first timestep after present

        self.always_valid = False  # Set to True if a single field should
                                   # be valid at all times

        self.is_lazy = False  # Generally False

        # Set projection for coordinate transformations
        self.simulation_SRS = False  # Avoid unnecessary vector rotation
        if self.proj is not None:
            self.projected = True
        else:
            if self.proj4 is not None:
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

    def set_convolution_kernel(self, convolve):
        """Set a convolution kernel or kernel size (of array of ones) used by `get_variables` on read variables."""
        self.convolve = convolve

    def prepare(self, extent, start_time, end_time):
        """Prepare reader for given simulation coverage in time and space."""
        logging.debug('Nothing to prepare for ' + self.name)
        pass  # to be overriden by specific readers

    def disabled_get_variables_interpolated(self, variables, profiles=None,
                                   profiles_depth=None, time=None,
                                   lon=None, lat=None, z=None,
                                   block=False, rotate_to_proj=None):
        """
        `get_variables_interpolated` is the interface to the basemodel, and is
        responsible for returning variables at the correct positions. This is done by:

            1. Calling `__get_variables__` which,
            2. calls `get_variables_derived`, which, finally
            3. calls `get_variables`.

        This function is responsible for setting up an interpolator,
        `ReaderBlock` for regular gridded datasets.

        Arguments:
            variables: string, or list of strings (standard_name) of
                requested variables. These must be provided by reader.

            profiles: N/A

            profiles_depth: N/A

            time: datetime or None, time at which data are requested.
                Can be None (default) if reader/variable has no time
                dimension (e.g. climatology or landmask).

            lon: N/A

            lat: N/A

            z: float or ndarray; vertical position (in meters, positive up)
                of requested points.
                default: 0 m (unless otherwise documented by reader)

            block: bool, see return below

            rotate_to_proj: N/A

          Returns:
            data: Dictionary
                keywords: variables (string)
                values:
                    - 1D ndarray of len(x) if block=False. Nearest values
                        (neighbour) of requested position are returned.
                    - 3D ndarray encompassing all requested points in
                        x,y,z domain if block=True. It is task of invoking
                        application (OpenDriftSimulation) to perform
                        interpolation in space and time.
        """

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
            env_before = self.get_variables_impl(variables, profiles,
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
                    self.get_variables_impl(blockvariables_before, profiles,
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
                        self.get_variables_impl(blockvariables_after, profiles,
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
            rx = np.float32(rx)
            ry = np.float32(ry)
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
            d = self.get_variables_interpolated(
                lon=np.atleast_1d(lon[i]), lat=np.atleast_1d(lat[i]), z=np.atleast_1d(0),
                time=closest_time, variables=variables, rotate_to_proj='+proj=latlong')[0]
            for var in variables:
                data[var][i] = d[var][0]

        return(data)
