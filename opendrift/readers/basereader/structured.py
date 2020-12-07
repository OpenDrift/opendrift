import numpy as np
from abc import abstractmethod

from opendrift.readers.interpolation.structured import ReaderBlock
from .variables import Variables
import logging
logger = logging.getLogger(__name__)

class StructuredReader(Variables):
    """
    A structured reader. Data is gridded on a regular grid. Used by e.g.:
    :class:`opendrift.readers.reader_netCDF_CF_generic.Reader`.

    .. seealso::

        :py:mod:`opendrift.readers`
    """
    var_block_before = None
    var_block_after  = None
    interpolation = 'linearNDFast'

    def __init__(self):
        super().__init__()

        # Dictionaries to store blocks of data for reuse (buffering)
        self.var_block_before = {}  # Data for last timestep before present
        self.var_block_after = {}   # Data for first timestep after present

    def _get_variables_interpolated_(self, variables, profiles,
                                   profiles_depth, time,
                                   reader_x, reader_y, z):
        # Find reader time_before/time_after
        time_nearest, time_before, time_after, i1, i2, i3 = \
            self.nearest_time(time)
        logger.debug('Reader time:\n\t\t%s (before)\n\t\t%s (after)' %
                      (time_before, time_after))

        # For variables which are not time dependent, we do not care about time
        static_variables = [
                'sea_floor_depth_below_sea_level',
                'land_binary_mask'
                ]
        if time == time_before or all(v in static_variables for v in variables):
            time_after = None

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
            reader_data_dict = \
                self._get_variables_impl_(blockvariables_before, profiles,
                                    profiles_depth, time_before,
                                    reader_x, reader_y, z)
            self.var_block_before[blockvars_before] = \
                ReaderBlock(reader_data_dict,
                            interpolation_horizontal=self.interpolation)
            try:
                len_z = len(self.var_block_before[blockvars_before].z)
            except:
                len_z = 1
            logger.debug(('Fetched env-block (size %ix%ix%i) ' +
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
                reader_data_dict = \
                    self._get_variables_impl_(blockvariables_after, profiles,
                                        profiles_depth, time_after,
                                        reader_x, reader_y, z)
                self.var_block_after[blockvars_after] = \
                    ReaderBlock(
                        reader_data_dict,
                        interpolation_horizontal=self.interpolation)
                try:
                    len_z = len(self.var_block_after[blockvars_after].z)
                except:
                    len_z = 1

                logger.debug(('Fetched env-block (size %ix%ix%i) ' +
                                'for time after (%s)') %
                                (len(self.var_block_after[blockvars_after].x),
                                len(self.var_block_after[blockvars_after].y),
                                len_z, time_after))
                block_after = self.var_block_after[blockvars_after]

        if (block_before is not None and block_before.covers_positions(
            reader_x, reader_y) is False) or (\
            block_after is not None and block_after.covers_positions(
                reader_x, reader_y) is False):
            logger.warning('Data block from %s not large enough to '
                            'cover element positions within timestep. '
                            'Buffer size (%s) must be increased.' %
                            (self.name, str(self.buffer)))

        ############################################################
        # Interpolate before/after blocks onto particles in space
        ############################################################
        self.timer_start('interpolation')
        logger.debug('Interpolating before (%s) in space  (%s)' %
                        (block_before.time, self.interpolation))
        env_before, env_profiles_before = block_before.interpolate(
                reader_x, reader_y, z, variables,
                profiles, profiles_depth)

        if (time_after is not None) and (time_before != time):
            logger.debug('Interpolating after (%s) in space  (%s)' %
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
        if (time_after is not None) and (time_before != time):
            weight_after = ((time - time_before).total_seconds() /
                            (time_after - time_before).total_seconds())
            logger.debug(('Interpolating before (%s, weight %.2f) and'
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
                logger.debug('Interpolating profiles in time')
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
            logger.debug('No time interpolation needed - right on time.')
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

        return env, env_profiles

