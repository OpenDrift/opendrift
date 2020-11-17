from .variables import Variables

class AnalyticalReader(Variables):
    """
    An analytical and continuous  (in space and time) reader, able to provide
    exact values at any desired point (within bounds).
    """
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

        # Analytical reader, continous in space and time
        self.timer_end('preparing')
        env = self._get_variables(variables, profiles,
                                            profiles_depth,
                                            time,
                                            #time_before,
                                            reader_x, reader_y, z,
                                            block=block)
        self.logger.debug('Fetched env-before')
        self.timer_start('preparing')

        if profiles is not None:
            if 'env_profiles_before' in locals():
                env_profiles = env_profiles_before
            else:
                # Copying data from environment to vertical profiles
                env_profiles = {'z': profiles_depth}
                for var in profiles:
                    env_profiles[var] = np.ma.array([env[var], env[var]])

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

