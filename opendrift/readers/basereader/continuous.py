import numpy as np
from .variables import Variables
from .consts import vector_pairs_xy

class ContinuousReader(Variables):
    """
    A continuous (in space and time) reader, able to provide
    exact values at any desired point (within bounds). This reader type is
    suitable for constant readers, analytical readers, or readers that are
    static and continuous within the valid domain (e.g. the landmask reader).
    """
    def _get_variables_interpolated_(self, variables, profiles,
                                   profiles_depth, time,
                                   lon, lat, z,
                                   block, rotate_to_proj,
                                   ind_covered, reader_x, reader_y):

        env = self.get_variables_impl(variables, profiles,
                                            profiles_depth,
                                            time,
                                            #time_before,
                                            reader_x, reader_y, z,
                                            block=block)

        self.logger.debug('Fetched env-before')
        env_profiles = None
        if profiles is not None:
            # Copying data from environment to vertical profiles
            env_profiles = {'z': profiles_depth}
            for var in profiles:
                env_profiles[var] = np.ma.array([env[var], env[var]])

        return env, env_profiles

