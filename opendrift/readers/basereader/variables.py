from datetime import timedelta
import numpy as np
import pyproj
from bisect import bisect_left
from abc import abstractmethod
import logging
logger = logging.getLogger(__name__)

from opendrift.timer import Timeable
from opendrift.errors import OutsideSpatialCoverageError, OutsideTemporalCoverageError, VariableNotCoveredError
from .consts import standard_names, vector_pairs_xy


class ReaderDomain(Timeable):
    """
    Projection, spatial and temporal domain of reader.
    """
    name = None

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
    shape = None

    ## Temporal
    start_time = None
    end_time = None
    time_step = None
    times = None

    """Setting this to `True` overrides temporal and spatial bounds checks.
    Also useful for readers that are constant and do not have a temporal
    dimension."""
    always_valid = False

    def center(self):
        """
        Returns center of reader (in lon, lat)
        """

        if any(xx is None for xx in [self.xmin, self.xmax, self.ymin, self.ymax]):
            return None, None

        x = self.xmin + (self.xmax - self.xmin) / 2
        y = self.ymin + (self.ymax - self.ymin) / 2
        lo, la = self.xy2lonlat(x, y)
        return(lo[0], la[0])

    def rotate_vectors(self, reader_x, reader_y, u_component, v_component,
                       proj_from, proj_to):
        if isinstance(u_component, list):  # Looping recursively over ensemble members
            uout = []
            vout = []
            for ucomp, vcomp in zip(u_component, v_component):
                ucomprot, vcomprot = self.rotate_vectors(
                    reader_x, reader_y, ucomp, vcomp, proj_from, proj_to)
                uout.append(ucomprot)
                vout.append(vcomprot)
            return uout, vout
        """Rotate vectors from one crs to another."""

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
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        if self.proj.crs.is_geographic:
            if 'ob_tran' in str(self.proj4):
                logger.debug('NB: Converting degrees to radians ' +
                             'due to ob_tran srs')
                x = np.radians(np.array(x))
                y = np.radians(np.array(y))
                return self.proj(x, y, inverse=True)
            else:
                return x, y
        else:
            return self.proj(x, y, inverse=True)

    def lonlat2xy(self, lon, lat):
        """
        Calculate lon,lat from given x,y (scalars/arrays) in own projection.
        """
        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)
        if 'ob_tran' in str(self.proj4):
            x, y = self.proj(lon, lat, inverse=False)
            return np.degrees(x), np.degrees(y)
        elif self.proj.crs.is_geographic:
            return lon, lat
        else:
            x, y = self.proj(lon, lat, inverse=False)
            return x, y

    def y_azimuth(self, lon, lat):
        """Calculate azimuth orientation of the y-axis of the reader CRS."""
        x0, y0 = self.lonlat2xy(lon, lat)
        distance = 1000.0  # Create points 1 km away to determine azimuth
        lon_2, lat_2 = self.xy2lonlat(x0, y0 + distance)
        geod = pyproj.Geod(ellps='WGS84')  # Define an ellipsoid
        y_az = geod.inv(lon_2, lat_2, lon, lat, radians=False)
        return y_az[0]

    def pixel_size(self):
        # Find typical pixel size (e.g. for calculating size of buffer)
        if self.delta_x is not None:
            pixelsize = self.delta_x
            if self.proj.crs.is_geographic is True or \
                    ('ob_tran' in self.proj4) or \
                    ('longlat' in self.proj4) or \
                    ('latlon' in self.proj4):
                pixelsize = pixelsize * 111000  # deg to meters
        else:
            pixelsize = None  # Pixel size not defined
        return pixelsize

    def _coverage_unit_(self):
        return "degrees"

    def __repr__(self):
        """String representation of the current reader."""
        outStr = '===========================\n'
        outStr += 'Reader: ' + self.name + '\n'
        outStr += 'Projection: \n  ' + self.proj4 + '\n'
        outStr += 'Coverage: [%s]\n' % self._coverage_unit_()
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
        if hasattr(self, 'realizations') and self.realizations is not None:
            outStr += 'Variables (%i ensemble members):\n' % len(
                self.realizations)
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

    def covers_positions_xy(self, x, y, z=0):
        """
        Return indices of input points covered by reader.

        Arguments in native projection of reader.
        """
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        z = np.atleast_1d(z)

        if self.global_coverage():
            # We need only check north-south and z coverage
            indices = np.where((y >= self.ymin) & (y <= self.ymax)
                               & (z >= self.zmin) & (z <= self.zmax))[0]
        else:
            xx = x

            if self.proj.crs.is_geographic:
                xx = self.modulate_longitude(x)

            indices = np.where((xx >= self.xmin) & (xx <= self.xmax)
                                & (y >= self.ymin) & (y <= self.ymax)
                                & (z >= self.zmin) & (z <= self.zmax))[0]

        try:
            return indices, x[indices], y[indices]
        except Exception as ex:
            logger.exception(ex)
            return indices, x, y

    def modulate_longitude(self, lons):
        """
        Modulate the input longitude to the domain supported by the reader.
        """

        yy = np.array([self.ymin, self.ymax, self.ymax, self.ymin])
        xx = np.array([self.xmin, self.xmin, self.xmax, self.xmax])

        exlons, _lats = self.xy2lonlat(xx, yy)

        if np.min(exlons) < 0:
            # Reader domain is from -180 to 180 or somewhere in between.
            assert np.max(exlons) <= 180.1

            lons = np.mod(lons+180, 360) - 180
        else:
            # Reader domain is from 0 to 360 or somewhere in between.
            assert np.max(exlons) <= 360.1

            lons = np.mod(lons, 360)

        return lons

    def covers_positions(self, lon, lat, z=0):
        """Return indices of input points covered by reader."""

        x, y = self.lonlat2xy(lon, lat)

        return self.covers_positions_xy(x, y, z)

    def global_coverage(self):
        """Return True if global coverage east-west"""

        dx = self.delta_x if self.delta_x is not None else 0
        if self.proj.crs.is_geographic is True:
            if (self.xmin - 2*dx <= 0) and (self.xmax + 2*dx >=
                                                    360):
                return True  # Global 0 to 360
            if (self.xmin - 2*dx <= -180) and (self.xmax + 2*dx
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
                raise VariableNotCoveredError('Variable not available: ' + variable +
                                 '\nAvailable parameters are: ' +
                                 str(self.variables))
        if (self.start_time is not None
                and time < self.start_time) and self.always_valid is False:
            raise OutsideTemporalCoverageError('Requested time (%s) is before first available '
                             'time (%s) of %s' % (time, self.start_time, self.name))
        if (self.end_time is not None
                and time > self.end_time) and self.always_valid is False:
            raise OutsideTemporalCoverageError('Requested time (%s) is after last available '
                             'time (%s) of %s' % (time, self.end_time, self.name))
        if self.global_coverage():
            outside = np.where(~np.isfinite(x + y) | (y < self.ymin)
                               | (y > self.ymax))[0]
        else:
            outside = np.where(~np.isfinite(x + y) | (x < self.xmin)
                               | (x > self.xmax) | (y < self.ymin)
                               | (y > self.ymax))[0]
        if np.size(outside) == np.size(x):
            lon, lat = self.xy2lonlat(x, y)
            raise OutsideSpatialCoverageError(('Argcheck: all %s particles (%.2f-%.2fE, ' +
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
            if len(self.times) > 1 and self.times[indx_before + 1] == time:
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

################################################################
# Methods to derive environment variables from others available
################################################################

def land_binary_mask_from_ocean_depth(env,in_name=None,out_name=None):
    env['land_binary_mask'] = np.float32(env['sea_floor_depth_below_sea_level'] <= 0)

def wind_from_speed_and_direction(env, in_name, out_name):
    if 'wind_from_direction' in env:
        wfd = env['wind_from_direction']
    else:
        wfd = -env['wind_to_direction']
    north_wind = -env['wind_speed']*np.cos(np.radians(wfd))
    east_wind = -env['wind_speed']*np.sin(np.radians(wfd))
    env['x_wind'] = east_wind
    env['y_wind'] = north_wind
    # Rotating might be necessary generally
    #x,y = np.meshgrid(env['x'], env['y'])
    #env['x_wind'], env['y_wind'] = self.rotate_vectors(
    #    x, y,
    #    east_wind, north_wind,
    #    None, self.proj)

def vector_from_speed_and_direction(env, in_name, out_name):
    # TODO: assuming here lonlat or Mercator projection  !!NB!!
    wfd = env[in_name[1]]
    env[out_name[0]] = env[in_name[0]]*np.cos(np.radians(wfd))
    env[out_name[1]] = env[in_name[0]]*np.sin(np.radians(wfd))

def reverse_direction(env, in_name, out_name):
    env[out_name[0]] = env[in_name[0]]

def magnitude_from_components(env, in_name, out_name):
    env[out_name[0]] = np.sqrt(
        env[in_name[0]]**2 + env[in_name[1]]**2)

################################################################

class Variables(ReaderDomain):
    """
    Handles reading and interpolation of variables.
    """
    variables = None
    derived_variables = None
    name = None

    buffer = 0

    def __init__(self):
        if self.derived_variables is None:
            self.derived_variables = {}
        if self.variables is None:
            self.variables = []

        # Deriving environment variables from other available variables
        self.environment_mappings = {
            #'wind_from_speed_and_direction': {
            #    'input': ['wind_speed', 'wind_from_direction'],
            #    'output': ['x_wind', 'y_wind'],
            #    'method': wind_from_speed_and_direction,
            #    #lambda reader, env: reader.wind_from_speed_and_direction(env)},
            #    'active': True},
            #'wind_from_speed_and_direction_to': {
            #    'input': ['wind_speed', 'wind_to_direction'],
            #    'output': ['x_wind', 'y_wind'],
            #    'method': wind_from_speed_and_direction,
            #    #lambda reader, env: reader.wind_from_speed_and_direction(env)},
            #    'active': True},
            'land_binary_mask_from_ocean_depth': {
                'input': ['sea_floor_depth_below_sea_level'],
                'output': ['land_binary_mask'],
                'method': land_binary_mask_from_ocean_depth,
                'active': False}
            }

        # Add automatically mappings from xcomp,ycomp <-> magnitude,direction
        for vector_pair in vector_pairs_xy:
            #if len(vector_pair) > 4:
                # TODO: temporarily disabling this test, as one test is failing,
                # as direction_to is automatically calculated from direction_from
                #self.environment_mappings[vector_pair[3]] = {
                #    'input': [vector_pair[4]],
                #    'output': [vector_pair[3]],
                #    'method': reverse_direction,
                #    'active': True
                #    }
                #self.environment_mappings[vector_pair[4]] = {
                #    'input': [vector_pair[3]],
                #    'output': [vector_pair[4]],
                #    'method': reverse_direction,
                #    'active': True
                #    }
            if len(vector_pair) >= 4:
                self.environment_mappings[str(vector_pair[0:2])] = {
                    'input': [vector_pair[2], vector_pair[3]],
                    'output': [vector_pair[0], vector_pair[1]],
                    'method': vector_from_speed_and_direction,
                    'active': True
                    }
            if len(vector_pair) > 2:
                self.environment_mappings[vector_pair[2]] = {
                    'input': [vector_pair[0], vector_pair[1]],
                    'output': [vector_pair[2]],
                    'method': magnitude_from_components,
                    'active': True
                    }

        super().__init__()

    def prepare(self, extent, start_time, end_time, max_speed):
        """Prepare reader for given simulation coverage in time and space."""
        logger.debug('Nothing more to prepare for ' + self.name)
        pass  # to be overriden by specific readers


    def activate_environment_mapping(self, mapping_name):
        if mapping_name not in self.environment_mappings:
            raise ValueError('Available environment mappings: ' + str(self.environment_mappings))

        em = self.environment_mappings[mapping_name]
        em['active'] = True
        if not all(item in em['output'] for item in self.variables) and \
                    all(item in self.variables for item in em['input']):
            for v in em['output']:
                if v not in self.variables:
                    logger.debug('Adding variable mapping: %s -> %s' % (em['input'], v))
                    self.variables.append(v)
                    self.derived_variables[v] = em['input']

    def __calculate_derived_environment_variables__(self, env):
        for m in self.environment_mappings:
            em = self.environment_mappings[m]
            if em['active'] is False:
                continue
            if not all(item in em['output'] for item in self.variables) and \
                    all(item in env for item in em['input']):
                for v in em['output']:
                    logger.debug('Calculating variable mapping: %s -> %s' % (em['input'], v))
                    em['method'](env, em['input'], em['output'])

    def set_buffer_size(self, max_speed, time_coverage=None):
        '''
        Adjust buffer to minimise data block size needed to cover elements.

        The buffer size is calculated from the maximum anticpated speed.
        Seeding over a large area or over longer time can easily cause
        particles to be located outside the block. This is not critical, but
        causes interpolation to be one-sided in time for the relevant
        particles.

        Args:

            max_speed (m/s): the maximum speed anticipated for particles.
            time_coverage (timedelta): the time span to cover

        '''
        self.buffer = 0
        pixelsize = self.pixel_size()
        if pixelsize is not None:
            if self.time_step is not None:
                time_step_seconds = self.time_step.total_seconds()
            else:
                if time_coverage is None:
                    logger.warning('Assuming time step of 1 hour for ' + self.name)
                    time_step_seconds = 3600  # 1 hour if not given
                else:
                    time_step_seconds = time_coverage.total_seconds()
            time_step_seconds = abs(time_step_seconds)
            self.buffer = int(
                np.ceil(max_speed * time_step_seconds / pixelsize)) + 2
            logger.debug('Setting buffer size %i for reader %s, assuming '
                         'a maximum average speed of %g m/s and time span of %s' %
                         (self.buffer, self.name, max_speed,timedelta(seconds=time_step_seconds)))

    def __check_env_coordinates__(self, env):
        """
        Make sure x and y are floats (and not e.g. int64)
        """
        if 'x' in env.keys():
            env['x'] = np.array(env['x'], dtype=np.float32)
            env['y'] = np.array(env['y'], dtype=np.float32)

    @staticmethod
    def __check_variable_array__(name, variable):
        """
        Ensure arrays are not masked arrays and that values are within valid ranges.
        """

        # Convert any masked arrays to NumPy arrays
        if isinstance(variable, np.ma.MaskedArray):
            variable = variable.astype(np.float32)
            variable = variable.filled(np.nan)

        # Mask values outside valid_min, valid_max (self.standard_names)
        if name in standard_names.keys():
            logger.debug("Checking %s for invalid values" % name)
            if isinstance(variable, list):
                logger.debug(
                    'Min-max checking for ensemble data (%i members)' %
                    (len(variable)))
                # Recursive
                return [
                    Variables.__check_variable_array__(name, v)
                    for v in variable
                ]
            # with np.errstate(invalid='ignore'):
            invalid_indices = np.logical_and(
                np.isfinite(variable),
                np.logical_or(variable < standard_names[name]['valid_min'],
                              variable > standard_names[name]['valid_max']))
            if np.sum(invalid_indices) > 0:
                invalid_values = variable[invalid_indices]
                logger.warning(
                    'Invalid values (%s to %s) found for %s, replacing with NaN'
                    % (invalid_values.min(), invalid_values.max(), name))
                logger.warning('(allowed range: [%s, %s])' %
                               (standard_names[name]['valid_min'],
                                standard_names[name]['valid_max']))
                variable[invalid_indices] = np.nan

        return variable

    def __check_env_arrays__(self, env):
        """
        Ensure arrays are not masked arrays and that values are within valid
        ranges.

        .. seealso::

            Disabled in `StructuredReader` because variables are valided before
            entered into interpolator:

            :meth:`.structured.StructuredReader.__check_env_arrays__`
        """
        variables = [
            var for var in env.keys() if var not in ['x', 'y', 'time']
        ]

        for variable in variables:
            env[variable] = self.__check_variable_array__(
                variable, env[variable])

        return env

    @abstractmethod
    def _get_variables_interpolated_(self, variables, profiles, profiles_depth,
                                     time, reader_x, reader_y, z):
        """
        This method _must_ be implemented by every reader. Usually by
        subclassing one of the reader types (e.g.
        :class:`structured.StructuredReader`).

        Arguments are in _native projection_ of reader.

        .. seealso:

            * :meth:`get_variables_interpolated_xy`.
            * :meth:`get_variables_interpolated`.
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

        .. seealso::

            * :meth:`get_variables_interpolated`.
            * :meth:`_get_variables_interpolated_`.
        """
        self.timer_start('total')
        # Raise error if time not not within coverage of reader
        if not self.covers_time(time):
            raise OutsideTemporalCoverageError('%s is outside time coverage (%s - %s) of %s' %
                             (time, self.start_time, self.end_time, self.name))

        self.timer_start('preparing')

        x = np.atleast_1d(x)
        numx = len(x)  # To check later if all points were covered

        y = np.atleast_1d(y)
        z = np.atleast_1d(z) if z is not None else np.zeros((1, ))

        assert numx == len(y), "x, y, and z must have the same length"
        assert len(z) == 1 or numx == len(
            z), "x, y, and z must have the same length"

        ind_covered, xx, yy = self.covers_positions_xy(x, y, z)
        if len(ind_covered) == 0:
            lon, lat = self.xy2lonlat(x, y)
            raise OutsideSpatialCoverageError(
                ('All %s particles (%.2f-%.2fE, %.2f-%.2fN) ' +
                 'are outside domain of %s (%s)') %
                 (len(x), lon.min(), lon.max(), lat.min(),
                 lat.max(), self.name, self.coverage_string()))
        x = xx
        y = yy

        self.timer_end('preparing')

        # Make copy of z to avoid modifying original array
        if len(z) == 1 and len(x) > 1:
            z = z.copy() * np.ones(x.shape)
        z = z.copy()[ind_covered]

        logger.debug(f'Fetching variables from {self.name} covering {len(ind_covered)} elements')
        self.timer_start('reading')

        # Filter derived variables
        derived = []
        derived_input = []

        for var in variables:
            if var in self.derived_variables:
                logger.debug("Scheduling %s to be derived from %s" %
                             (var, self.derived_variables[var]))
                derived_input.extend(self.derived_variables[var])
                derived.append(var)

        for v in derived:
            variables.remove(v)
        variables.extend(list(set(derived_input)))

        env, env_profiles = self._get_variables_interpolated_(
            variables, profiles, profiles_depth, time, x, y, z)

        # Calculate derived variables
        if len(derived) > 0:
            self.__calculate_derived_environment_variables__(env)
            variables.extend(derived)

        for e in [env, env_profiles]:
            if e is not None:
                self.__check_env_coordinates__(e)
                self.__check_env_arrays__(e)

        self.timer_end('reading')

        # Rotating vectors fields
        if rotate_to_proj is not None:
            if self.proj.crs.is_geographic and 'ob_tran' not in self.proj4:
                logger.debug('Reader projection is latlon - '
                             'rotation of vectors is not needed.')
            else:
                vector_pairs = []
                for var in variables:
                    for vector_pair in vector_pairs_xy:
                        if var in vector_pair[0:2]:
                            counterpart = list(set(vector_pair[0:2]) -
                                               set([var]))[0]
                            if counterpart in variables:
                                vector_pairs.append(vector_pair[0:2])
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
        `get_variables_interpolated` is the main interface to
        :class:`opendrift.basemodel.OpenDriftSimulation`, and is responsible
        for returning variables at the correct positions.

        Readers should implement :meth:`_get_variables_interpolated_`.

        Arguments:
            variables: string, or list of strings (standard_name) of
                requested variables. These must be provided by reader.

            profiles: List of variable names that should be returned for the range in `profiles_depth`.

            profiles_depth: Profiles variables will be retrieved from surface and down to this depth. The exact z-depth are given by the reader and returned as `z` variable in `env_profiles`.

            time: datetime or None, time at which data are requested.
                Can be None (default) if reader/variable has no time
                dimension (e.g. climatology or landmask).

            lon: longitude, 1d array.

            lat: latitude, 1d array, same length as lon.

            z: float or ndarray; vertical position (in meters, positive up)
                of requested points. either scalar or same length as lon, lat.
                default: 0 m (unless otherwise documented by reader)

            block: bool, see return below

            rotate_to_proj: N/A

          Returns:

            (env, env_profiles)

            Interpolated variables at x, y and z. `env` contains values at a fixed depth (`z`), while `env_profiles` contains depth-profiles in the range `profile_depth` for the variables listed in `profiles` for each element (in `x`, `y`). The exact depth is determined by the reader and specified in
            `env_profiles['z']`. Thus variables in `env_profiles` are not interpolated in z-direction.

        .. seealso::

            :meth:`get_variables_interpolated_xy`.

        """
        assert set(variables).issubset(self.variables), f"{variables} is not subset of {self.variables}"

        lon = self.modulate_longitude(lon)
        x, y = self.lonlat2xy(lon, lat)

        env, env_profiles = self.get_variables_interpolated_xy(
            variables, profiles, profiles_depth, time, x, y, z, rotate_to_proj)

        return env, env_profiles
