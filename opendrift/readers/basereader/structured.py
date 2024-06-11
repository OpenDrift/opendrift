import pickle
import numpy as np
import pyproj
from scipy.ndimage import map_coordinates
from abc import abstractmethod
from pathlib import Path

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

    # TODO: should the variables below not be instance variables, and not class variables?
    clipped = 0
    x = None
    y = None
    interpolation = 'linearNDFast'
    convolve = None  # Convolution kernel or kernel size
    # set these in reader to save interpolators to file
    save_interpolator = None
    interpolator_filename = None
    

    # Used to enable and track status of parallel coordinate transformations.
    __lonlat2xy_parallel__ = None
    __disable_parallel__ = False

    def __init__(self):
        if self.proj is None and (self.proj4 is None
                                  or self.proj4 == 'fakeproj'):

            logger.warning(
                "No proj string or projection could be derived, using 'fakeproj'. This assumes that the variables are structured and gridded approximately equidistantly on the surface (i.e. in meters). This must be guaranteed by the user. You can get rid of this warning by supplying a valid projection to the reader."
            )

            from scipy.interpolate import LinearNDInterpolator
            import copy
            from . import fakeproj

            # `projected` is set to True if `fakeproj` is used
            self.projected = None
            self.shape = None

            self.proj4 = 'None'
            self.proj = fakeproj.fakeproj()
            self.projected = False
            
            self.xmin = self.ymin = 0.
            self.delta_x = self.delta_y = 1.
            self.xmax = self.lon.shape[1] - 1
            self.ymax = self.lon.shape[0] - 1
            self.numx = self.xmax
            self.numy = self.ymax
            self.x = np.arange(0, self.xmax+1)
            self.y = np.arange(0, self.ymax+1)
    
            # Making interpolator (lon, lat) -> x
            # save to speed up next time
            if self.save_interpolator and self.interpolator_filename is not None:
                interpolator_filename = Path(self.interpolator_filename).with_suffix('.pickle')
            else:
                interpolator_filename = f'{self.name}_interpolators.pickle'
            
            if self.save_interpolator and Path(interpolator_filename).is_file():
                logger.info('Loading previously saved interpolator for lon,lat to x,y conversion.')
                with open(interpolator_filename, 'rb') as file_handle:
                    interp_dict = pickle.load(file_handle)
                    spl_x = interp_dict["spl_x"]
                    spl_y = interp_dict["spl_y"]
            
            else:
                logger.info('Making interpolator for lon,lat to x,y conversion...')

                block_x, block_y = np.mgrid[self.xmin:self.xmax + 1,
                                            self.ymin:self.ymax + 1]
                block_x, block_y = block_x.T, block_y.T
            
                spl_x = LinearNDInterpolator(
                    (self.lon.ravel(), self.lat.ravel()),
                    block_x.ravel(),
                    fill_value=np.nan)
                # Reusing x-interpolator (deepcopy) with data for y
                spl_y = copy.deepcopy(spl_x)
                spl_y.values[:, 0] = block_y.ravel()
                # Call interpolator to avoid threading-problem:
                # https://github.com/scipy/scipy/issues/8856
                spl_x((0, 0)), spl_y((0, 0))
                
                if self.save_interpolator:
                    logger.info('Saving interpolator for lon,lat to x,y conversion.')

                    interp_dict = {"spl_x": spl_x, "spl_y": spl_y}
                    with open(interpolator_filename, 'wb') as f:
                        pickle.dump(interp_dict, f)
            
            self.spl_x = spl_x
            self.spl_y = spl_y


        else:
            self.projected = True

        super().__init__()

        # Dictionaries to store blocks of data for reuse (buffering)
        self.var_block_before = {}  # Data for last timestep before present
        self.var_block_after = {}  # Data for first timestep after present

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

    def prepare(self, extent, start_time, end_time, max_speed):
        """Prepare reader for given simulation coverage in time and space."""
        logger.debug('Clearing cache for reader %s before starting new simulation' % self.name)
        self.var_block_before = {}
        self.var_block_after = {}
        if self.time_step is None and start_time is not None:
            # Set buffer large enough for whole simulation
            logger.debug('Time step is None for %s, setting buffer size large nough for whole simulation' % self.name)
            self.set_buffer_size(max_speed, end_time-start_time)
        else:
            self.set_buffer_size(max_speed, self.time_step)

        super().prepare(extent, start_time, end_time, max_speed)

    def set_convolution_kernel(self, convolve):
        """Set a convolution kernel or kernel size (of array of ones) used by `get_variables` on read variables."""
        self.convolve = convolve

    def __convolve_block__(self, env):
        """
        Convolve arrays with a kernel, if reader.convolve is set
        """
        if self.convolve is not None:
            from scipy import ndimage
            N = self.convolve
            if isinstance(N, (int, np.integer)):
                kernel = np.ones((N, N))
                kernel = kernel / kernel.sum()
            else:
                kernel = N
            logger.debug('Convolving variables with kernel: %s' % kernel)
            for variable in env:
                if variable in ['x', 'y', 'z', 'time']:
                    pass
                else:
                    if env[variable].ndim == 2:
                        env[variable] = ndimage.convolve(env[variable],
                                                         kernel,
                                                         mode='nearest')
                    elif env[variable].ndim == 3:
                        env[variable] = ndimage.convolve(env[variable],
                                                         kernel[:, :, None],
                                                         mode='nearest')
        return env

    def lon_range(self):
        if not self.global_coverage():
            raise ValueError('Only valid for readers with global coverage')
        if self.xmin < 0:
            return '-180to180'
        else:
            return '0to360'

    def _get_variables_interpolated_(self, variables, profiles, profiles_depth,
                                     time, reader_x, reader_y, z):

        # For global readers, we shift coordinates to match actual lon range
        if self.global_coverage():
            if self.lon_range() == '-180to180':
                logger.debug('Shifting coordinates to -180-180')
                reader_x = np.mod(reader_x + 180, 360) - 180
            elif self.lon_range() == '0to360':
                logger.debug('Shifting coordinates to 0-360')
                reader_x = np.mod(reader_x, 360)
        elif self.proj.crs.is_geographic and self.xmin>0:
            logger.debug('Modulating longitudes to 0-360 for self.name')
            reader_x = np.mod(reader_x, 360)

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
            mz = np.append(z, [0, -profiles_depth])
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
                    self.__convolve_block__(
                self.get_variables(blockvariables_before, time_before,
                                    mx, my, mz)
                )
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
                self.var_block_after[blockvars_after] = block_before
            else:
                reader_data_dict = self.__convolve_block__(
                    self.get_variables(blockvariables_after, time_after, mx,
                                       my, mz))
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
                           'Buffer size (%s) must be increased. See `Variables.set_buffer_size`.' %
                           (self.name, str(self.buffer)))
            # TODO: could add dynamic increase of buffer size here

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
        if (time_after is not None) and (time_before != time) and self.always_valid is False:
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
                env_profiles['z'] = env_profiles_before['z'][0:numlayers]
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
                    env_profiles = {'z': [0, -profiles_depth]}
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
            y = np.atleast_1d(y)
            x = np.atleast_1d(x)

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
            self.__lonlat2xy_parallel__ = False
            return super().lonlat2xy(lon, lat)
        else:
            # For larger arrays, we split and calculate in parallel
            num_elements = len(np.atleast_1d(lon))
            if num_elements > 10000 and not self.__disable_parallel__:
                from multiprocessing import cpu_count
                from concurrent.futures import ThreadPoolExecutor

                self.__lonlat2xy_parallel__ = True

                nproc = cpu_count()
                logger.debug('Running lonlat2xy in parallel using %d threads' %
                             nproc)

                # Chunk arrays
                split_lon = np.array_split(lon, nproc)
                split_lat = np.array_split(lat, nproc)

                with ThreadPoolExecutor() as x:
                    out_x = np.concatenate(
                        list(x.map(self.spl_x, zip(split_lon, split_lat))))
                    out_y = np.concatenate(
                        list(x.map(self.spl_y, zip(split_lon, split_lat))))

                return (out_x, out_y)

            else:
                logger.debug('Calculating lonlat2xy sequentially')
                self.__lonlat2xy_parallel__ = False
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

    def get_ocean_depth_area_volume(self, lonmin, lonmax, latmin, latmax):
        """Get depth, area and volume of ocean basin within given coordinates"""

        # Extract ocean depth within given boundaries
        background = 'sea_floor_depth_below_sea_level'
        rx, ry = self.lonlat2xy([lonmin, lonmax, lonmax, lonmin], [latmin, latmin, latmax, latmax])
        rx = np.linspace(rx.min(), rx.max(), 10)
        ry = np.linspace(ry.min(), ry.max(), 10)
        data = self.get_variables(background, time=None, x=rx, y=ry)

        x, y = np.meshgrid(data['x'], data['y'])
        lon, lat = self.xy2lonlat(x, y)

        depth = data[background]
        depth = np.ma.masked_where(lon<lonmin, depth)
        depth = np.ma.masked_where(lon>lonmax, depth)
        depth = np.ma.masked_where(lat<latmin, depth)
        depth = np.ma.masked_where(lat>latmax, depth)

        volume = np.nansum(depth*self.pixel_size()*self.pixel_size())
        area = volume/np.nanmean(depth)

        return np.nanmin(depth), np.nanmax(depth), np.nanmean(depth), area, volume

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

    def _slice_variable_(self,
                         var,
                         indxTime=None,
                         indy=None,
                         indx=None,
                         indz=None,
                         indrealization=None):
        """
        Slice variable depending on number of dimensions available.

        Args:

            All arguments can be `slice` objects or index.

        Returns:

            `var` sliced using the slices or indexes necessary to use depending
            on number of dimensions available.

        Raises:

            Unsupported number of dimensions (outside 2..5) raises an exception.

        """
        # NOTE: use match expressions when PEP-634 (Py 3.10) is (widely)
        #       available.
        if var.ndim == 2:
            return var[indy, indx]
        elif var.ndim == 3:
            return var[indxTime, indy, indx]
        elif var.ndim == 4:
            return var[indxTime, indz, indy, indx]
        elif var.ndim == 5:  # Ensemble data
            return var[indxTime, indz, indrealization, indy, indx]
        else:
            raise Exception('Wrong dimension of variable: %s: %d' %
                            (var, var.ndim))
