import sys
import importlib
import logging
from bisect import bisect_left
from abc import abstractmethod, ABCMeta
from scipy.interpolate import RectBivariateSpline, LinearNDInterpolator
from scipy.spatial import KDTree
import numpy as np


try:
    import pyproj  # Import pyproj
except:
    try:
        # ...alternatively use version included with Basemap
        from mpl_toolkits.basemap import pyproj
    except:
        raise ImportError('pyproj needed for coordinate transformations,'
                          ' please install from '
                          'https://code.google.com/p/pyproj/')

# Som valid (but extreme) ranges for checking that values are reasonable
standard_names = {
    'x_wind': {'valid_min': -50, 'valid_max': 50},
    'y_wind': {'valid_min': -50, 'valid_max': 50},
    'x_sea_water_velocity': {'valid_min': -10, 'valid_max': 10},
    'y_sea_water_velocity': {'valid_min': -10, 'valid_max': 10} }

# Identify x-y vector components/pairs for rotation (NB: not east-west pairs!)
vector_pairs_xy = [
    ['x_wind', 'y_wind'],
    ['x_sea_water_velocity', 'y_sea_water_velocity']]

class Reader(object):
    """Parent Reader class, to be subclassed by specific readers.
    """

    __metaclass__ = ABCMeta

    return_block = True  # By default, all readers should be
                         # cabable of returning blocks of data

    start_time = None
    # Dictionaries to store blocks of data for reuse (buffering)
    var_block_before = {}  # Data for last timestep before present
    var_block_after = {}   # Data for first timestep after present

    # Mapping variable names, e.g. from east-north to x-y, temporarily
    # presuming coordinate system then is lon-lat for equivalence
    variable_aliases = {'eastward_sea_water_velocity': 'x_sea_water_velocity',
                        'northward_sea_water_velocity': 'y_sea_water_velocity'}

    def __init__(self):
        # Common constructor for all readers

        # Set projection for coordinate transformations
        self.proj = pyproj.Proj(self.proj4)

        # Check if there are holes in time domain
        if self.start_time is not None:
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

    @abstractmethod
    def get_variables(self, variables, time=None,
                      x=None, y=None, depth=None, block=False):
        """Method which must be invoked by any reader (subclass).

        Obtain and return values of the requested variables at all positions
        (x, y, depth) closest to given time.
        
        Arguments:
            variables: string, or list of strings (standard_name) of
                requested variables. These must be provided by reader.
            time: datetime or None, time at which data are requested.
                Can be None (default) if reader/variable has no time 
                dimension (e.g. climatology or landmask).
            x, y: float or ndarrays; coordinates of requested points in the
                Spatial Reference System (SRS) of the reader (NB!!)
            depth: float or ndarray; depth (in meters) of requested points.
                default: 0 m (unless otherwise documented by reader)
            block: bool, see return below

          Returns:
            data: Dictionary
                keywords: variables (string)
                values:
                    - 1D ndarray of len(x) if block=False. Nearest values
                        (neichbour) of requested position are returned.
                    - 3D ndarray encompassing all requested points in
                        x,y,depth domain if block=True. It is task of invoking
                        application (OpenDriftSimulation) to perform 
                        interpolation in space and time.
        """

    def interpolate_block(self, block, x, y, depth=None):
        """Interpolating a 2D or 3D block onto given positions."""

        env = {}
        for var in block.keys():
            if not hasattr(block[var], 'ndim'):
                continue
            if block[var].ndim == 2:
                ## Create spline for interpolation
                block_x, block_y = np.meshgrid(block['x'], block['y'])
                data = block[var]
                data[data.mask] = np.nan
                spl = LinearNDInterpolator(
                        (block_y.ravel(), block_x.ravel()), data.ravel())
                env[var] = spl(y, x)  # Evaluate at positions

                ## Make and use KDTree for interpolation
                #tree = KDTree(zip(block['x'].ravel(), block['y'].ravel()))
                #d,p = tree.query(zip(x,y), k=1) #nearest point, gives distance and point index
                ### NB - or transpose array first ?
                #env[var] = block[var].ravel()[p]

            elif block[var].ndim == 3:
                raise ValueError('Not yet implemented')
        return env

    def get_variables_from_buffer(self, variables, time=None,
                                  lon=None, lat=None, depth=None,
                                  block=False, rotate_to_proj=None):
        """Wrapper around get_variables(), reading from buffer if available.
        
        Also performs interpolation in the following order:
        - horizontally (x-y, bilinear)
        - vertically (z/depth, TBD)
        - time (linear)
        
        """

        # Raise error if not not within time
        if not self.covers_time(time):
            raise ValueError('Outside time coverage of ' + self.name)
        # Find reader time_before/time_after 
        time_nearest, time_before, time_after, i1, i2, i3 = \
            self.nearest_time(time)
        logging.debug('Reader time:\n\t%s (nearest)\n\t%s '
                      '(before)\n\t%s (after)' %
                      (time_nearest, time_before, time_after))
        if time == time_before:
            time_after = None
        # Check which particles are covered (indep of time) 
        ind_covered = self.covers_positions(lon, lat, depth)
        if len(ind_covered)==0:
            raise ValueError('All particles are outside domain '
                             'of ' + self.name)  

        reader_x, reader_y = self.lonlat2xy(lon, lat)
        reader_x_min = reader_x.min()
        reader_x_max = reader_x.max()
        reader_y_min = reader_y.min()
        reader_y_max = reader_y.max()

        if block is False or self.return_block is False:
            env_before = self.get_variables(variables, time_before,
                            reader_x, reader_y, depth,
                            block=block)
            logging.debug('Fetched env-before')
            if time_after is not None:
                env_after = self.get_variables(variables, time_after,
                                reader_x, reader_y, depth,
                                block=block)
                logging.debug('Fetched env-after')
            else:
                env_after = None
        else:
            # Swap before- and after-blocks if matching times
            if self.var_block_before.has_key(str(variables)):
                block_before_time = self.var_block_before[
                                        str(variables)]['time']
                if self.var_block_after.has_key(str(variables)):
                    block_after_time = self.var_block_after[
                                        str(variables)]['time']
                    if block_before_time != time_before:
                        if block_after_time == time_before:
                            self.var_block_before[
                                str(variables)]  = \
                            self.var_block_after[
                                str(variables)]
                    if block_after_time != time_after:
                        if block_before_time == time_before:
                            self.var_block_after[
                                str(variables)]  = \
                            self.var_block_before[
                                str(variables)]  
            # Fetch data, if no buffer is available
            if not self.var_block_before.has_key(str(variables)) or \
                self.var_block_before[str(variables)]['time'] != time_before:
                self.var_block_before[str(variables)] = \
                    self.get_variables(variables, time_before,
                        reader_x, reader_y, depth,
                        block=block)
                logging.debug('Fetched env-block for time before (%s)' %
                        (time_before))
            if not self.var_block_after.has_key(str(variables)) or \
                self.var_block_after[str(variables)]['time'] != time_after:
                if time_after is None:
                    self.var_block_after[str(variables)] = \
                        self.var_block_before[str(variables)]
                else:
                    self.var_block_after[str(variables)] = \
                        self.get_variables(variables, time_after,
                            reader_x, reader_y, depth,
                            block=block)
                    logging.debug('Fetched env-block for time after (%s)' %
                        (time_after))
            # check if buffer-block covers these particles
            x_before = self.var_block_before[str(variables)]['x']
            y_before = self.var_block_before[str(variables)]['y']
            x_after = self.var_block_after[str(variables)]['x']
            y_after = self.var_block_after[str(variables)]['y']
            if (reader_x_min < x_before.min() or
                reader_x_max > x_before.max() or
                reader_y_min < y_before.min() or
                reader_y_max > y_before.max()):
                logging.debug('Some elements not covered by before-block')
                update  # to be implemented
            else:
                logging.debug('All elements covered by before-block')
            if (reader_x_min < x_after.min() or
                reader_x_max > x_after.max() or
                reader_y_min < y_after.min() or
                reader_y_max > y_after.max()):
                logging.debug('Some elements not covered by after-block')
                update  # to be implemented
            else:
                logging.debug('All elements covered by after-block')
            # Interpolate before/after onto particles in space
            #   x-y
            logging.debug('Interpolating before (%s) in space' %
                (self.var_block_before[str(variables)]['time']))
            env_before = self.interpolate_block(
                            self.var_block_before[str(variables)],
                            reader_x, reader_y, depth)
            logging.debug('Interpolating after (%s) in space' %
                (self.var_block_after[str(variables)]['time']))
            env_after = self.interpolate_block(
                            self.var_block_after[str(variables)],
                            reader_x, reader_y, depth)
            #   depth interpolation - TBD

        # Time interpolation
        if (time_after is not None) and (time_before != time):
            weight = ((time - time_before).total_seconds() /
                       (time_after - time_before).total_seconds())
            logging.debug('Interpolating before (%s, weight %f) and '
                          'after (%s, weight %f) in time' %
                          (self.var_block_before[str(variables)]['time'],
                          1 - weight,
                          self.var_block_after[str(variables)]['time'],
                          weight))
            env = {}
            for var in variables:
                # Weighting together, and masking invalid entries
                env[var] = np.ma.masked_invalid(
                            (env_before[var] * (1 - weight) + 
                            env_after[var] * weight))

                if var in standard_names.keys():
                    if (env[var].min() < standard_names[var]['valid_min']) \
                        or (env[var].max() > standard_names[var]['valid_max']):
                            logging.info('Invalid values found for ' + var)
                            logging.ingo(env[var])
                            sys.exit('quitting') 
        else:
            logging.debug('No time interpolation needed - right on time.')
            env = env_before

        # Rotate vectors
        if rotate_to_proj is not None:
            vector_pairs = []
            for var in variables:
                for vector_pair in vector_pairs_xy:
                    if var in vector_pair:
                        counterpart = list(set(vector_pair) - set([var]))[0]
                        if counterpart in variables:
                            vector_pairs.append(vector_pair)
                        else:
                            sys.exit('Missing component of vector pair:' +
                                        counterpart)
            # Extract unique vector pairs
            vector_pairs = [list(x) for x in set(tuple(x)
                                for x in vector_pairs)]
            # Calculate azimuth angle between coordinate systems (y-axis)
            if self.delta_y is None:
                delta_y = 1000
            else:
                delta_y = self.delta_y*.1  # 10% of grid separation
            x2, y2 = pyproj.transform(self.proj, rotate_to_proj,
                                        reader_x, reader_y)
            x2_delta, y2_delta = pyproj.transform(self.proj,
                                    rotate_to_proj,
                                    reader_x, reader_y + delta_y)
            rot_angle_rad = np.arctan2(x2_delta - x2, y2_delta - y2)
            logging.debug('Rotating vectors to srs: "%s"' %
                            rotate_to_proj.srs)
            for vector_pair in vector_pairs:
                logging.debug('Rotating %s an angle of %f to %f degrees' %
                    (str(vector_pair), np.degrees(rot_angle_rad).min(),
                        np.degrees(rot_angle_rad).min()))
                u2 = (env[vector_pair[0]]*np.cos(rot_angle_rad) -
                        env[vector_pair[1]]*np.sin(rot_angle_rad))
                v2 = (env[vector_pair[0]]*np.sin(rot_angle_rad) +
                        env[vector_pair[1]]*np.cos(rot_angle_rad))
                env[vector_pair[0]] = u2
                env[vector_pair[1]] = v2

        return env

    def xy2lonlat(self, x, y):
        """Calculate x,y in own projection from given lon,lat (scalars/arrays).
        """
        if self.proj.is_latlong():
            return x,y
        else:
            return self.proj(x, y, inverse=True)

    def lonlat2xy(self, lon, lat):
        """Calculate lon,lat from given x,y (scalars/arrays) in own projection.
        """
        if self.proj.is_latlong():
            return lon,lat
        else:
            return self.proj(lon, lat, inverse=False)

    def y_azimuth(self, lon, lat):
        """Calculate the azimuth orientation of the y-axis of the reader SRS."""
        x0, y0 = self.lonlat2xy(lon, lat)
        distance = 1000.0  # Create points 1 km away to determine azimuth
        lon_2, lat_2 = self.xy2lonlat(x0, y0 + distance)
        geod = pyproj.Geod(ellps='WGS84')  # Define an ellipsoid
        y_az = geod.inv(lon_2, lat_2, lon, lat, radians=False)
        dist = y_az[2]
        return y_az[0]

    def covers_time(self, time):

        if self.start_time is None:
            return True  # No time limitations of reader
        if (time < self.start_time) or (time > self.end_time):
            return False
        else:
            return True

    def covers_positions(self, lon, lat, depth):
        """Return indices of input points covered by reader."""

        # Compensate for wrapping about 0 or 180 longitude
        if self.proj.is_latlong():
            if self.xmax > 180:
                lon[lon<0] = lon[lon<0] + 360

        # Calculate x,y coordinates from lon,lat
        x, y = self.lonlat2xy(lon, lat)

        indices = np.where((x > self.xmin) & (x < self.xmax) &
                           (y > self.xmin) & (y < self.ymax))[0]

        return indices

    def check_coverage(self, time, lon, lat):
        """Check which points are within coverage of reader.
        
        Checks that requested positions and time are within coverage of
        this reader, and that it can provide the requested variable(s).
        Returns the input arguments, possibly modified/corrected (below)

        Arguments:
            See function get_variables for definition.

        Returns:
            x, y: coordinates of point in spatial reference system of reader
            indices: indices of the input points which are inside domain

        Raises:
            ValueError:
                - if requested time is outside coverage of reader.
                - if all requested positions are outside coverage of reader.
        """

        # Check time
        if self.start_time is not None and (time is not None and
                                            time < self.start_time):
            raise ValueError('Requested time (%s) is before first available '
                             'time (%s) of %s' % (time, self.start_time,
                                                  self.name))
        if self.end_time is not None and (time is not None and
                                          time > self.end_time):
            raise ValueError('Requested time (%s) is after last available '
                             'time (%s) of %s' % (time, self.end_time,
                                                  self.name))

        # Compensate for wrapping about 0 or 180 longitude
        if self.proj.is_latlong():
            if self.xmax > 180:
                lon[lon<0] = lon[lon<0] + 360

        # Calculate x,y coordinates from lon,lat
        x, y = self.lonlat2xy(lon, lat)

        indices = np.where((x > self.xmin) & (x < self.xmax) &
                           (y > self.xmin) & (y < self.ymax))[0]
        if len(indices)==0:
            raise ValueError('All particles are outside domain '
                             'of ' + self.name)

        return x[indices], y[indices], indices

    def check_arguments(self, variables, time, x, y, depth):
        """Check validity of arguments input to method get_variables.
        
        Checks that requested positions and time are within coverage of
        this reader, and that it can provide the requested variable(s).
        Returns the input arguments, possibly modified/corrected (below)

        Arguments:
            See function get_variables for definition.

        Returns:
            variables: same as input, but converted to list if given as string.
            time: same as input, or start_time of reader if given as None.
            x, y, depth: same as input, but converted to ndarrays
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
            time = self.start_time  # Get data from first timestep, if not given

        # Convert variables to list and x,y to ndarrays
        if isinstance(variables, str):
            variables = [variables]
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        depth = np.asarray(depth)

        for variable in variables:
            if variable not in self.variables:
                raise ValueError('Variable not available: ' + variable +
                                 '\nAvailable parameters are: ' +
                                 str(self.variables))
        if self.start_time is not None and time < self.start_time:
            raise ValueError('Requested time (%s) is before first available '
                             'time (%s) of %s' % (time, self.start_time,
                                                  self.name))
        if self.end_time is not None and time > self.end_time:
            raise ValueError('Requested time (%s) is after last available '
                             'time (%s) of %s' % (time, self.end_time,
                                                  self.name))
        outside = np.where((x < self.xmin) | (x > self.xmax ) |
                           (y < self.xmin) | (y > self.ymax ))
        if np.size(outside) == np.size(x):
            raise ValueError('All particles are outside domain '
                             'of ' + self.name)

        return variables, time, x, y, depth, outside

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

        if self.start_time is None:
            return None, None, None, None, None, None
        if hasattr(self, 'times'):  # Time as array, possibly with holes
            indx_before = bisect_left(self.times, time) - 1
            time_before = self.times[indx_before]
            indx_after = indx_before + 1
            time_after = self.times[indx_after]
            if (time - time_before) < (time_after - time):
                indx_nearest = indx_before
            else:
                indx_nearest = indx_after
            nearest_time = self.times[indx_nearest]
        else:  # Time step is constant (no holes)
            indx = float((time - self.start_time).total_seconds()) / \
                float(self.time_step.total_seconds())
            indx_nearest = int(round(indx))
            nearest_time = self.times[indx_nearest]
            indx_before = int(np.floor(indx))
            time_before = self.times[indx_before]
            indx_after = int(np.ceil(indx))
            time_after = self.times[indx_after]
        return nearest_time, time_before, time_after,\
                indx_nearest, indx_before, indx_after

    def index_of_closest_depths(self, requestedDepths):
        """Return (internal) index of depths closest to requested depths.

        Depths of layers (of ocean model) are not assumed to be constant.
        """
        ind_depth = [np.abs(np.subtract.outer(
            self.depths, requestedDepths)).argmin(0)]
        return ind_depth, self.depths[ind_depth]

    def indices_min_max_depth(self, depths):
        """
        Return min and max indices of internal depth dimension,
        covering the requested depths. Needed when block is requested (True).

        Arguments:
            depths: ndarray of floats, in meters
        """
        minIndex = (self.depths <= depths.min()).argmin() - 1
        maxIndex = (self.depths >= depths.max()).argmax()
        return minIndex, maxIndex

    def __repr__(self):
        """String representation of the current reader."""
        outStr = '===========================\n'
        outStr += 'Reader: ' + self.name + '\n'
        outStr += 'Projection: \n  ' + self.proj4 + '\n'
        outStr += 'Coverage: \n'
        outStr += '  xmin: %f   xmax: %f   step: %f\n' % \
            (self.xmin, self.xmax, self.delta_x or 0)
        outStr += '  ymin: %f   ymax: %f   step: %f\n' % \
            (self.ymin, self.ymax, self.delta_y or 0)
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
        if hasattr(self, 'depths'):
            np.set_printoptions(suppress=True)
            outStr += 'Depths [m]: \n  ' + str(self.depths) + '\n'
        outStr += 'Available time range:\n'
        outStr += '  start: ' + str(self.start_time) + \
                  '   end: ' + str(self.end_time) + \
                  '   step: ' + str(self.time_step) + '\n'
        if self.start_time is not None:
            outStr += '    %i times (%i missing)\n' % (
                        self.expected_time_steps, self.missing_time_steps)
        outStr += 'Variables:\n'
        for variable in self.variables:
            outStr += '  ' + variable + '\n'
        outStr += '===========================\n'
        return outStr

    def plot(self):
        """Plot geographical coverage of reader."""

        try:
            from mpl_toolkits.basemap import Basemap
            import matplotlib.pyplot as plt
            from matplotlib.patches import Polygon
        except:
            sys.exit('Basemap is needed to plot coverage map.')

        # Initialise map, stereographic projection centred on domain
        x0 = (self.xmin + self.xmax) / 2
        y0 = (self.ymin + self.ymax) / 2
        lon0, lat0 = self.xy2lonlat(x0, y0)
        width = np.max([self.xmax-self.xmin, self.ymax-self.ymin])*1.5
        map = Basemap(projection='stere', resolution='i',
                      lat_ts=lat0, lat_0=lat0, lon_0=lon0,
                      width=width, height=width)
        map.drawcoastlines()
        map.fillcontinents(color='coral')
        map.drawparallels(np.arange(-90.,90.,5.))
        map.drawmeridians(np.arange(-180.,181.,5.))
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
        mapproj = pyproj.Proj(map.proj4string)
        xm, ym = mapproj(lon, lat)
        #map.plot(xm, ym, color='gray')
        boundary = Polygon(zip(xm, ym), alpha=0.5, ec='k', fc='b')
        plt.gca().add_patch(boundary)
# add patch to the map
        plt.title(self.name)
        plt.xlabel('Time coverage: %s to %s' % (self.start_time, self.end_time))
        plt.show()
