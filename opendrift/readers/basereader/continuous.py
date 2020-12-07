import numpy as np
from .variables import Variables
from .consts import vector_pairs_xy
import logging
logger = logging.getLogger(__name__)

class ContinuousReader(Variables):
    """
    A continuous (in space and time) reader, able to provide
    exact values at any desired point (within bounds). This reader type is
    suitable for constant readers, analytical readers, or readers that are
    static and continuous within the valid domain (e.g. the landmask reader).

    .. seealso::

        :py:mod:`opendrift.readers`
    """
    def _get_variables_interpolated_(self, variables, profiles,
                                   profiles_depth, time,
                                   reader_x, reader_y, z):

        env = self._get_variables_impl_(variables, profiles,
                                            profiles_depth,
                                            time,
                                            reader_x, reader_y, z)

        logger.debug('Fetched env-before')
        env_profiles = None
        if profiles is not None:
            # Copying data from environment to vertical profiles
            env_profiles = {'z': profiles_depth}
            for var in profiles:
                env_profiles[var] = np.ma.array([env[var], env[var]])

        return env, env_profiles

