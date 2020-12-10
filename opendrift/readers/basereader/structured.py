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

    Attributes:

        projected: is `True` if :class:`.fakeproj.fakeproj` is used because of missing projection information. The data points are assumed to be approximately equidistant on the surface (i.e. in meters).

        clipped: pixels to to remove along boundary (e.g. in case of bad data).

    .. seealso::

        :py:mod:`opendrift.readers`
    """
    # `projected` is set to True if `fakeproj` is used
    projected = None
    shape = None
    clipped = 0
    x = None
    y = None

    var_block_before = None
    var_block_after = None
    interpolation = 'linearNDFast'

    def __init__(self):
        if self.proj is None and (self.proj4 is None
                                  or self.proj4 == 'fakeproj'):

            logger.warning(
                "No proj string or projection could be derived, using 'fakeproj'. This assumes that the variables are structured and gridded approximately equidistantly on the surface (i.e. in meters). This must be guaranteed by the user. You can get rid of this warning by suppling a valid projection to the reader."
            )

            from scipy.interpolate import LinearNDInterpolator
            import copy
            from . import fakeproj

            self.proj4 = 'None'
            self.proj = fakeproj.fakeproj()
            self.projected = False
            logger.info('Making interpolator for lon,lat to x,y conversion...')
            self.xmin = self.ymin = 0.
            self.delta_x = self.delta_y = 1.
            self.xmax = self.lon.shape[1] - 1
            self.ymax = self.lon.shape[0] - 1
            self.numx = self.xmax
            self.numy = self.ymax

            block_x, block_y = np.mgrid[self.xmin:self.xmax + 1,
                                        self.ymin:self.ymax + 1]
            block_x, block_y = block_x.T, block_y.T

            # Making interpolator (lon, lat) -> x
            self.spl_x = LinearNDInterpolator(
                (self.lon.ravel(), self.lat.ravel()),
                block_x.ravel(),
                fill_value=np.nan)
            # Reusing x-interpolator (deepcopy) with data for y
            self.spl_y = copy.deepcopy(self.spl_x)
            self.spl_y.values[:, 0] = block_y.ravel()
            # Call interpolator to avoid threading-problem:
            # https://github.com/scipy/scipy/issues/8856
            self.spl_x((0, 0)), self.spl_y((0, 0))
        else:
            self.projected = True

        super().__init__()

        # Dictionaries to store blocks of data for reuse (buffering)
        self.var_block_before = {}  # Data for last timestep before present
        self.var_block_after = {}   # Data for first timestep after present

    @abstractmethod
    def get_variables(self, variables, time=None, x=None, y=None, z=None):
        """
        Obtain a _block_ of values of the requested variables at all positions
        (x, y, z) closest to given time. These will be stored in
        :class:`opendrift.readers.interpolation.structured.ReaderBlock` and
        accessed from there.

        Arguments:
            variables: list of variables.

            time: datetime or None, time at which data are requested.

            x, y: float or ndarrays; coordinates of requested points.

            z: float or ndarray; vertical position (in meters, positive up)

          Returns:
            Dictionary

            keywords: variables (string)
            values: 2D ndarray bounding x and y.
        """

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
            'sea_floor_depth_below_sea_level', 'land_binary_mask'
        ]
        if time == time_before or all(v in static_variables
                                      for v in variables):
            time_after = None


        if profiles is not None:
            # If profiles are requested for any parameters, we
            # add two fake points at the end of array to make sure that the
            # requested block has the depth range required for profiles
            mx = np.append(reader_x, [reader_x[-1], reader_x[-1]])
            my = np.append(reader_y, [reader_y[-1], reader_y[-1]])
            mz = np.append(z, [profiles_depth[0], profiles_depth[1]])
        else:
            mx = reader_x
            my = reader_y
            mz = z

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
                self.get_variables(blockvariables_before, time_before,
                                    mx, my, mz)
            self.var_block_before[blockvars_before] = \
                ReaderBlock(reader_data_dict,
                            interpolation_horizontal=self.interpolation)
            try:
                len_z = len(self.var_block_before[blockvars_before].z)
            except:
                len_z = 1
            logger.debug(
                ('Fetched env-block (size %ix%ix%i) ' + 'for time before (%s)')
                % (len(self.var_block_before[blockvars_before].x),
                   len(self.var_block_before[blockvars_before].y), len_z,
                   time_before))
            block_before = self.var_block_before[blockvars_before]
        if block_after is None or block_after.time != time_after:
            if time_after is None:
                self.var_block_after[blockvars_after] = \
                    block_before
            else:
                reader_data_dict = \
                    self.get_variables(blockvariables_after, time_after,
                                        mx, my, mz)
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
            reader_x, reader_y, z, variables, profiles, profiles_depth)

        if (time_after is not None) and (time_before != time):
            logger.debug('Interpolating after (%s) in space  (%s)' %
                         (block_after.time, self.interpolation))
            env_after, env_profiles_after = block_after.interpolate(
                reader_x, reader_y, z, variables, profiles, profiles_depth)

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
                env[var] = np.ma.masked_invalid(
                    (env_before[var] * (1 - weight_after) +
                     env_after[var] * weight_after))
            # Interpolating vertical profiles in time
            if profiles is not None:
                env_profiles = {}
                logger.debug('Interpolating profiles in time')
                # Truncating layers not present both before and after
                numlayers = np.minimum(len(env_profiles_before['z']),
                                       len(env_profiles_after['z']))
                env_profiles['z'] = env_profiles_before['z'][0:numlayers + 1]
                for var in env_profiles_before.keys():
                    if var == 'z':
                        continue
                    env_profiles_before[var] = np.atleast_2d(
                        env_profiles_before[var])
                    env_profiles_after[var] = np.atleast_2d(
                        env_profiles_after[var])
                    env_profiles[var] = (
                        env_profiles_before[var][0:numlayers, :] *
                        (1 - weight_after) +
                        env_profiles_after[var][0:numlayers, :] * weight_after)
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

    def __check_env_arrays__(self, env):
        """
        For the StructuredReader the variables are checked before entered into
        the ReaderBlock interpolator. This methods makes the second check a
        no-op.

        .. seealso::

            :meth:`.variables.Variables.__check_env_arrays__`.
        """
        return env

    def xy2lonlat(self, x, y):
        if self.projected:
            return super().xy2lonlat(x, y)
        else:
            np.seterr(invalid='ignore')  # Disable warnings for nan-values
            y = np.atleast_1d(np.array(y))
            x = np.atleast_1d(np.array(x))

            # NB: mask coordinates outside domain
            x[x < self.xmin] = np.nan
            x[x > self.xmax] = np.nan
            y[y < self.ymin] = np.nan
            y[y < self.ymin] = np.nan

            lon = map_coordinates(self.lon, [y, x],
                                  order=1,
                                  cval=np.nan,
                                  mode='nearest')
            lat = map_coordinates(self.lat, [y, x],
                                  order=1,
                                  cval=np.nan,
                                  mode='nearest')
            return (lon, lat)

    def lonlat2xy(self, lon, lat):
        if self.projected:
            return super().lonlat2xy(lon, lat)
        else:
            # Two methods needed for parallelisation
            def get_x(lon_part, lat_part, num):
                out_x[num] = self.spl_x(lon_part, lat_part)

            def get_y(lon_part, lat_part, num):
                out_y[num] = self.spl_y(lon_part, lat_part)

            # TODO: Use ThreadPool, maximum number of threads should not exceed
            # available cpus.

            # For larger arrays, we split and calculate in parallel
            # The number of CPUs to use can be improved/optimised
            num_elements = len(np.atleast_1d(lon))
            if num_elements > 100 and not hasattr(self,
                                                  'multiprocessing_fail'):
                try:
                    nproc = 2
                    if num_elements > 1000:
                        nproc = 16
                    if num_elements > 100000:
                        nproc = 32
                    if num_elements > 1000000:
                        nproc = 64
                    cpus = cpu_count()
                    nproc = np.minimum(nproc, cpus - 1)
                    nproc = np.maximum(2, nproc)
                    logger.debug(
                        'Running lonlat2xy in parallel, using %i of %i CPUs' %
                        (nproc, cpus))
                    split_lon = np.array_split(lon, nproc)
                    split_lat = np.array_split(lat, nproc)
                    out_x = Manager().dict()
                    out_y = Manager().dict()
                    processes = []
                    for i in range(nproc):
                        processes.append(
                            Process(target=get_x,
                                    args=(split_lon[i], split_lat[i], i)))
                        processes.append(
                            Process(target=get_y,
                                    args=(split_lon[i], split_lat[i], i)))
                    [p.start() for p in processes]
                    [p.join() for p in processes]
                    x = np.concatenate(out_x)
                    y = np.concatenate(out_y)
                    logger.debug('Completed lonlat2xy in parallel')
                    return (x, y)
                except Exception as e:
                    logger.warning('Parallelprocessing failed:')
                    logger.warning(e)
                    self.multiprocessing_fail = True
            else:
                if hasattr(self, 'multiprocessing_fail'):
                    logger.warning(
                        'Multiprocessing has previously failed, reverting to using single processor for lonlat -> xy'
                    )
                else:
                    # For smaller arrays, we run sequentially
                    pass
            logger.debug('Calculating lonlat->xy sequentially')
            x = self.spl_x(lon, lat)
            y = self.spl_y(lon, lat)
            return (x, y)

    def pixel_size(self):
        if self.projected:
            return super().pixel_size()
        else:
            lons, lats = self.xy2lonlat([self.xmin, self.xmax],
                                        [self.ymin, self.ymin])
            geod = pyproj.Geod(ellps='WGS84')  # Define an ellipsoid
            dist = geod.inv(lons[0], lats[0], lons[1], lats[1],
                            radians=False)[2]
            pixelsize = dist / self.shape[0]
            return pixelsize

    def _coverage_unit_(self):
        if self.projected:
            return super()._coverage_unit_()
        else:
            return "pixels"

    def _bbox_(self, x, y):
        """
        Find bounding box on grid containing points (x, y)
        """
        ix = (x - self.xmin) / self.delta_x
        ix0, ix1 = np.min(ix), np.max(ix)

        iy = (y - self.ymin) / self.delta_y
        iy0, iy1 = np.min(iy), np.max(iy)

        ix0 = np.max((self.clipped, ix0 - self.buffer)).astype(int)
        iy0 = np.max((self.clipped, iy0 - self.buffer)).astype(int)

        ix1 = np.min((self.numx - self.clipped, ix1 + self.buffer)).astype(int)
        iy1 = np.min((self.numy - self.clipped, iy1 + self.buffer)).astype(int)

        return (ix0, ix1, iy0, iy1)

    def _make_projected_grid_(self, lon, lat, eq_eps=1.e-1):
        """
        Make the projected grid in cases where `lon` and `lat` are present as
        2D variables, but not `x` and `y` and assert that it is approximately
        equidistant.

        Args:

            eq_eps: tolerance for equidistance checks.
        """

        if self.x is not None or self.y is not None:
            logger.error("x and y variables already exist!")

        logger.debug("Finding bounds of reader")
        assert len(lon.shape) == 2
        assert len(lat.shape) == 2

        self.X, self.Y = self.lonlat2xy(lon, lat)
        self.xmin, self.xmax = np.min(self.X[:]), np.max(self.X[:])
        self.ymin, self.ymax = np.min(self.Y[:]), np.max(self.Y[:])

        self.delta_x = np.diff(self.X).flat[0]
        self.delta_y = np.diff(self.Y, axis=0).flat[0]

        self.x = self.X[0, :]
        self.y = self.Y[:, 0]
        self.numx = len(self.x)
        self.numy = len(self.y)

        self.__validate_projected_grid__(eq_eps)

    def __validate_projected_grid__(self, eq_eps=1.e-1):
        """
        Validate that the projected grid is approximately equidistant.

        Args:

            eq_eps: tolerance for equidistance checks.

        Raises:

            AssertionError if not equidistant within `eq_eps`.
        """

        assert np.all(np.abs(self.delta_x - np.diff(self.X)) < eq_eps
                      ), "Grid is not equidistant in X direction"
        assert np.all(np.abs(self.delta_y - np.diff(self.Y, axis=0)) < eq_eps
                      ), "Grid is not equidistant in Y direction"

        assert np.all(
            np.abs(np.tile(self.x, (self.X.shape[0], 1)) - self.X) < eq_eps
        ), "X coordinates are not aligned along Y direction"
        assert np.all(
            np.abs(
                np.tile(np.atleast_2d(self.y).T, (1, self.Y.shape[1])) - self.Y
            ) < eq_eps), "Y coordinates are not aligned along X direction"
