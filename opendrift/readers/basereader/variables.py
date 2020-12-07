import numpy as np
import pyproj
from scipy.ndimage import map_coordinates
from bisect import bisect_left
from multiprocessing import Process, Manager, cpu_count
import copy
from abc import abstractmethod
import logging
logger = logging.getLogger(__name__)

from opendrift.timer import Timeable
from .consts import standard_names, vector_pairs_xy


class ReaderDomain(Timeable):
    """
    Projection, spatial and temporal domain of reader.
    """
    simulation_SRS = False
    projected = None
    proj4 = None
    proj = None

    lon = None
    lat = None
    xmin = None
    xmax = None
    ymin = None
    ymax = None
    zmin = -np.inf
    zmax = np.inf
    delta_x = None
    delta_y = None

    ## Temporal
    start_time = None
    end_time = None
    time_step = None
    times = None

    def rotate_vectors(self, reader_x, reader_y, u_component, v_component,
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
        x2_delta, y2_delta = transformer.transform(reader_x,
                                                   reader_y + delta_y)

        if proj_to.crs.is_geographic:
            geod = pyproj.Geod(ellps='WGS84')
            rot_angle_vectors_rad = np.radians(
                geod.inv(x2, y2, x2_delta, y2_delta)[0])
        else:
            rot_angle_vectors_rad = np.arctan2(x2_delta - x2, y2_delta - y2)
        logger.debug('Rotating vectors between %s and %s degrees.' %
                     (np.degrees(rot_angle_vectors_rad).min(),
                      np.degrees(rot_angle_vectors_rad).max()))
        rot_angle_rad = -rot_angle_vectors_rad
        u_rot = (u_component * np.cos(rot_angle_rad) -
                 v_component * np.sin(rot_angle_rad))
        v_rot = (u_component * np.sin(rot_angle_rad) +
                 v_component * np.cos(rot_angle_rad))

        return u_rot, v_rot

    def xy2lonlat(self, x, y):
        """Calculate x,y in own projection from given lon,lat (scalars/arrays).
        """
        if self.projected is True:
            if self.proj.crs.is_geographic:
                if 'ob_tran' in self.proj4:
                    logger.debug('NB: Converting degrees to radians ' +
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

    def y_azimuth(self, lon, lat):
        """Calculate azimuth orientation of the y-axis of the reader SRS."""
        x0, y0 = self.lonlat2xy(lon, lat)
        distance = 1000.0  # Create points 1 km away to determine azimuth
        lon_2, lat_2 = self.xy2lonlat(x0, y0 + distance)
        geod = pyproj.Geod(ellps='WGS84')  # Define an ellipsoid
        y_az = geod.inv(lon_2, lat_2, lon, lat, radians=False)
        dist = y_az[2]
        return y_az[0]

    def pixel_size(self):
        # Find typical pixel size (e.g. for calculating size of buffer)
        if self.projected is True:
            if self.delta_x is not None:
                pixelsize = self.delta_x
                if self.proj.crs.is_geographic is True or \
                        ('ob_tran' in self.proj4) or \
                        ('longlat' in self.proj4) or \
                        ('latlon' in self.proj4):
                    pixelsize = pixelsize * 111000  # deg to meters
            else:
                pixelsize = None  # Pixel size not defined
        else:
            lons, lats = self.xy2lonlat([self.xmin, self.xmax],
                                        [self.ymin, self.ymin])
            typicalsize = lons
            geod = pyproj.Geod(ellps='WGS84')  # Define an ellipsoid
            dist = geod.inv(lons[0], lats[0], lons[1], lats[1],
                            radians=False)[2]
            pixelsize = dist / self.shape[0]
        return pixelsize

    def covers_positions_xy(self, x, y, z=0):
        """
        Return indices of input points covered by reader.

        Arguments in native projection of reader.
        """
        z = np.atleast_1d(z)

        if self.global_coverage():
            # We need only check north-south and z coverage
            indices = np.where((y >= self.ymin) & (y <= self.ymax)
                               & (z >= self.zmin) & (z <= self.zmax))[0]
        else:
            indices = np.where((x >= self.xmin) & (x <= self.xmax)
                               & (y >= self.ymin) & (y <= self.ymax)
                               & (z >= self.zmin) & (z <= self.zmax))[0]

        try:
            return indices, x[indices], y[indices]
        except Exception as ex:
            logger.exception(ex)
            return indices, x, y

    def covers_positions(self, lon, lat, z=0):
        """Return indices of input points covered by reader."""

        x, y = self.lonlat2xy(lon, lat)

        return self.covers_positions_xy(x, y, z)

    def global_coverage(self):
        """Return True if global coverage east-west"""

        if self.proj.crs.is_geographic is True and self.delta_x is not None:
            if (self.xmin - self.delta_x <= 0) and (self.xmax + self.delta_x >=
                                                    360):
                return True  # Global 0 to 360
            if (self.xmin - self.delta_x <= -180) and (self.xmax + self.delta_x
                                                       >= 180):
                return True  # Global -180 to 180

        return False

    def domain_grid(self, npoints=1000):
        """Return arrays of lon,lat points spread evenly over reader domain."""
        numx = np.floor(np.sqrt(npoints))
        numy = np.floor(np.sqrt(npoints))
        x = np.linspace(self.xmin, self.xmax - self.delta_x, numx)
        y = np.linspace(self.ymin, self.ymax - self.delta_y, numy)
        x, y = np.meshgrid(x, y)
        lons, lats = self.xy2lonlat(x, y)
        return lons, lats

    def coverage_string(self):
        """Coverage of reader to be reported as string for debug output"""
        corners = self.xy2lonlat([self.xmin, self.xmin, self.xmax, self.xmax],
                                 [self.ymax, self.ymin, self.ymax, self.ymin])
        return '%.2f-%.2fE, %.2f-%.2fN' % (np.min(
            corners[0]), np.max(corners[0]), np.min(
                corners[1]), np.max(corners[1]))

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
        if isinstance(variables, str):
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
        if (self.start_time is not None
                and time < self.start_time) and self.always_valid is False:
            raise ValueError('Requested time (%s) is before first available '
                             'time (%s) of %s' %
                             (time, self.start_time, self.name))
        if (self.end_time is not None
                and time > self.end_time) and self.always_valid is False:
            raise ValueError('Requested time (%s) is after last available '
                             'time (%s) of %s' %
                             (time, self.end_time, self.name))
        if self.global_coverage():
            outside = np.where(~np.isfinite(x + y) | (y < self.ymin)
                               | (y > self.ymax))[0]
        else:
            outside = np.where(~np.isfinite(x + y) | (x < self.xmin)
                               | (x > self.xmax) | (y < self.ymin)
                               | (y > self.ymax))[0]
        if np.size(outside) == np.size(x):
            lon, lat = self.xy2lonlat(x, y)
            raise ValueError(('Argcheck: all %s particles (%.2f-%.2fE, ' +
                              '%.2f-%.2fN) are outside domain of %s (%s)') %
                             (len(lon), lon.min(), lon.max(), lat.min(),
                              lat.max(), self.name, self.coverage_string()))

        return variables, time, x, y, z, outside

    def covers_time(self, time):
        if self.always_valid is True:
            return True
        if self.start_time is None:
            return True  # No time limitations of reader
        if (time < self.start_time) or (time > self.end_time):
            return False
        else:
            return True

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
        if self.times is not None:  # Time as array, possibly with holes
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
            nearest_time = self.start_time + indx_nearest * self.time_step
            indx_before = int(np.floor(indx))
            time_before = self.start_time + indx_before * self.time_step
            indx_after = int(np.ceil(indx))
            time_after = self.start_time + indx_after * self.time_step
        return nearest_time, time_before, time_after,\
            indx_nearest, indx_before, indx_after


class Variables(ReaderDomain):
    """
    Handles reading and interpolation of variables.
    """
    derived_variables = None
    name = None

    buffer = 0
    convolve = None  # Convolution kernel or kernel size

    environment_mappers = []
    environment_mappings = {
        'wind_from_speed_and_direction': {
            'input': ['wind_speed', 'wind_to_direction'],
            'output': ['x_wind', 'y_wind'],
            'method':
            lambda reader, var: reader.wind_from_speed_and_direction(var)
        },
        'testvar': {
            'input': ['sea_ice_thickness'],
            'output': ['istjukkleik']
        }
    }

    def __init__(self):
        self.derived_variables = {}

    def set_convolution_kernel(self, convolve):
        """Set a convolution kernel or kernel size (of array of ones) used by `get_variables` on read variables."""
        self.convolve = convolve

    def set_buffer_size(self, max_speed, max_vertical_speed=None):
        '''Adjust buffer to minimise data block size needed to cover elements'''
        self.buffer = 0
        pixelsize = self.pixel_size()
        if pixelsize is not None:
            if self.time_step is not None:
                time_step_seconds = self.time_step.total_seconds()
            else:
                time_step_seconds = 3600  # 1 hour if not given
            self.buffer = np.int(
                np.ceil(max_speed * time_step_seconds / pixelsize)) + 2
            logger.debug('Setting buffer size %i for reader %s, assuming '
                         'a maximum average speed of %g m/s.' %
                         (self.buffer, self.name, max_speed))

    def wind_from_speed_and_direction(self, env):
        north_wind = env['wind_speed'] * np.cos(
            np.radians(env['wind_to_direction']))
        east_wind = env['wind_speed'] * np.sin(
            np.radians(env['wind_to_direction']))
        env['x_wind'] = east_wind
        env['y_wind'] = north_wind
        # Rotating might be necessary generally
        #x,y = np.meshgrid(env['x'], env['y'])
        #env['x_wind'], env['y_wind'] = self.rotate_vectors(
        #    x, y,
        #    east_wind, north_wind,
        #    None, self.proj)

    def calculate_derived_environment_variables(self, env):

        if 'x_wind' in self.derived_variables and 'wind_speed' in env.keys():
            self.wind_from_speed_and_direction(env)

    def _get_variables_impl_(self, variables, profiles, profiles_depth, time, x,
                           y, z):
        """Wrapper around reader-specific function get_variables()

        Performs some common operations which should not be duplicated:

            - monitor time spent by this reader
            - convert any numpy arrays to masked arrays
            - add points to x, y and z to get profiles if necessary
            - convolve vectors if `self.convolve` is specified

        This function calls :meth:`__get_variables_derived__` which eventually
        calls :meth:`get_variables` on the implementing reader.
        """
        logger.debug('Fetching variables from ' + self.name)
        self.timer_start('reading')

        if profiles is not None:
            # If profiles are requested for any parameters, we
            # add two fake points at the end of array to make sure that the
            # requested block has the depth range required for profiles
            x = np.append(x, [x[-1], x[-1]])
            y = np.append(y, [y[-1], y[-1]])
            z = np.append(z, [profiles_depth[0], profiles_depth[1]])

        env = self.__get_variables_derived__(variables, time, x, y, z)

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
                    logger.warning(
                        'Skipping min-max checking for ensemble data')
                    continue
                with np.errstate(invalid='ignore'):
                    invalid_indices = np.logical_and(
                        np.isfinite(env[variable]),
                        np.logical_or(
                            env[variable] <
                            standard_names[variable]['valid_min'],
                            env[variable] >
                            standard_names[variable]['valid_max']))
                if np.sum(invalid_indices) > 0:
                    invalid_values = env[variable][invalid_indices]
                    logger.warning(
                        'Invalid values (%s to %s) found for %s, replacing with NaN'
                        %
                        (invalid_values.min(), invalid_values.max(), variable))
                    logger.warning('(allowed range: [%s, %s])' %
                                   (standard_names[variable]['valid_min'],
                                    standard_names[variable]['valid_max']))
                    env[variable][invalid_indices] = np.nan

        # Convolve arrays with a kernel, if reader.convolve is set
        if self.convolve is not None:
            from scipy import ndimage
            N = self.convolve
            if isinstance(N, (int, np.integer)):
                kernel = np.ones((N, N))
                kernel = kernel / kernel.sum()
            else:
                kernel = N
            logger.debug('Convolving variables with kernel: %s' % kernel)
            for variable in env.keys():
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

        self.timer_end('reading')

        return env

    def __get_variables_derived__(self, variables, *args, **kwargs):
        """Wrapper around get_variables, adding derived"""
        if isinstance(variables, str):
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

    @abstractmethod
    def get_variables(self, variables, time=None, x=None, y=None, z=None):
        """
        Method which must be implemented by all reader-subclasses.

        .. warning::

            Warning: In the future this method is likely to be a requirement of
            each reader, and it will be up to the reader-type how values are
            retrieved.

        Obtain and return values of the requested variables at all positions
        (x, y, z) closest to given time.

        Arguments:
            variables: string, or list of strings (standard_name) of
            requested variables. These must be provided by reader.

            time: datetime or None, time at which data are requested.
            Can be None (default) if reader/variable has no time
            dimension (e.g. climatology or landmask).

            x, y: float or ndarrays; coordinates of requested points in the
            Spatial Reference System (SRS) _of the reader (NB!!)_.

            z: float or ndarray; vertical position (in meters, positive up)
            of requested points. default: 0 m (unless otherwise documented by reader)

            block: bool, see return below

          Returns:
            data: Dictionary

            keywords: variables (string)

            values:

            - 1D ndarray of len(x) if StructuredReader. Nearest values
                (neighbour) of requested position are returned.

            - 3D ndarray encompassing all requested points in
                x,y,z domain if UnstructuredReader. It is task of invoking
                application (OpenDriftSimulation) to perform
                interpolation in space and time.
        """

    @abstractmethod
    def _get_variables_interpolated_(self, variables, profiles, profiles_depth,
                                     time, reader_x, reader_y, z):
        """
        Implemented by different reader types (e.g. :class:`structured.StructuredReader`).

        Arguments are in native projection of reader.
        """
        pass

    def get_variables_interpolated_xy(self,
                                      variables,
                                      profiles=None,
                                      profiles_depth=None,
                                      time=None,
                                      x=None,
                                      y=None,
                                      z=None,
                                      rotate_to_proj=None):
        """
        Get variables in native projection of reader.
        """
        self.timer_start('total')
        # Raise error if time not not within coverage of reader
        if not self.covers_time(time):
            raise ValueError('%s is outside time coverage (%s - %s) of %s' %
                             (time, self.start_time, self.end_time, self.name))

        self.timer_start('preparing')

        x = np.atleast_1d(x)
        y = np.atleast_1d(y)

        numx = len(x)  # To check later if all points were covered

        ind_covered, xx, yy = self.covers_positions_xy(x, y, z)
        if len(ind_covered) == 0:
            logger.error("All particles outside domain!")
            raise ValueError(('All %s particles (%.2f-%.2fE, %.2f-%.2fN) ' +
                              'are outside domain of %s (%s)') %
                             (len(x), x.min(), x.max(), y.min(), y.max(),
                              self.name, self.coverage_string()))
        x = xx
        y = yy

        self.timer_end('preparing')

        # Make copy of z to avoid modifying original array
        if len(z) == 1 and len(x) > 1:
            z = z.copy() * np.ones(x.shape)
        z = z.copy()[ind_covered]

        env, env_profiles = self._get_variables_interpolated_(
            variables, profiles, profiles_depth, time, x, y, z)

        # Rotating vectors fields
        if rotate_to_proj is not None:
            if self.simulation_SRS is True:
                logger.debug('Reader SRS is the same as calculation SRS - '
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
                                logger.warning(
                                    'Missing component of vector pair, cannot rotate:'
                                    + counterpart)
                # Extract unique vector pairs
                vector_pairs = [
                    list(x) for x in set(tuple(x) for x in vector_pairs)
                ]

                if len(vector_pairs) > 0:
                    self.timer_start('rotating vectors')
                    for vector_pair in vector_pairs:
                        env[vector_pair[0]], env[vector_pair[1]] = \
                            self.rotate_vectors(x, y,
                                                env[vector_pair[0]],
                                                env[vector_pair[1]],
                                                self.proj, rotate_to_proj)
                        if profiles is not None and vector_pair[0] in profiles:
                            env_profiles[vector_pair[0]], env_profiles[vector_pair[1]] = \
                                    self.rotate_vectors (x, y,
                                            env_profiles[vector_pair[0]],
                                            env_profiles[vector_pair[1]],
                                            self.proj,
                                            rotate_to_proj)

                    self.timer_end('rotating vectors')

        # Masking non-covered pixels
        self.timer_start('masking')
        if len(ind_covered) != numx:
            logger.debug('Masking %i elements outside coverage' %
                         (numx - len(ind_covered)))
            for var in variables:
                tmp = np.nan * np.ones(numx)
                tmp[ind_covered] = env[var].copy()
                env[var] = np.ma.masked_invalid(tmp)
                # Filling also in missing columns
                # for env_profiles outside coverage
                if env_profiles is not None and var in env_profiles.keys():
                    tmp = np.nan * np.ones((env_profiles[var].shape[0], numx))
                    tmp[:, ind_covered] = env_profiles[var].copy()
                    env_profiles[var] = np.ma.masked_invalid(tmp)

        self.timer_end('masking')
        self.timer_end('total')

        return env, env_profiles

    def get_variables_interpolated(self,
                                   variables,
                                   profiles=None,
                                   profiles_depth=None,
                                   time=None,
                                   lon=None,
                                   lat=None,
                                   z=None,
                                   rotate_to_proj=None):
        """
        `get_variables_interpolated` is the interface to
        :class:`opendrift.basemodel.OpenDriftSimulation`, and is responsible
        for returning variables at the correct positions. This is done by:

            1. Calling :meth:`_get_variables_interpolated_` which,
            2. calls :meth:`_get_variables_impl_`, which
            3. calls :meth:`__get_variables_derived__`, which
            4. calls :meth:`get_variables`.

        :meth:`_get_variables_impl_`: Works on every variable. If profiles_depth, adds a point at start and end in order to get a full block. This seems specific to `StructuredReader`. Needs to work on both env and env_profiles, but also modifies the behavior to make env_profiles work in the first place.

        :meth:`__get_variables_derived__`: Calculates derived variables from variables present in reader. Needs to work on both env and env_profiles.

        Arguments:
            variables: string, or list of strings (standard_name) of
                requested variables. These must be provided by reader.

            profiles: List of variable names that should be returned for the range in `profiles_depth`.

            profiles_depth: A range [z-start, z-end] for which to return values for profile-variables. The exact z-depth are given by the reader and returned as `z` variable in `env_profiles`.

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

            (env, env_profiles)

            Interpolated variables at x, y and z. `env` contains values at a fixed depth (`z`), while `env_profiles` contains depth-profiles in the range `profile_depth` for the variables listed in `profiles` for each element (in `x`, `y`). The exact depth is determined by the reader and specified in
            `env_profiles['z']`. Thus variables in `env_profiles` are not interpolated in z-direction.

        """

        x, y = self.lonlat2xy(lon, lat)

        env, env_profiles = self.get_variables_interpolated_xy(
            variables, profiles, profiles_depth, time, x, y, z, rotate_to_proj)

        return env, env_profiles
