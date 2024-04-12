from abc import abstractmethod
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

    @abstractmethod
    def get_variables(self, variables, time=None, x=None, y=None, z=None):
        """
        Obtain and return values of the requested variables at all positions
        (x, y, z) closest to given time.

        Returns:

            Dictionary with arrays of length len(x) with values at exact positions x, y and z.
        """

    def _get_variables_interpolated_(self, variables, profiles,
                                   profiles_depth, time,
                                   reader_x, reader_y, z):

        env = self.get_variables(variables, time, reader_x, reader_y, z)

        logger.debug('Fetched env-before')
        env_profiles = None
        if profiles is not None:
            # Copying data from environment to vertical profiles
            env_profiles = {'z': [0, -profiles_depth]}
            for var in profiles:
                env_profiles[var] = np.ma.array([env[var], env[var]])

        return env, env_profiles

