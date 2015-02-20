import importlib
from abc import abstractmethod, ABCMeta
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


class Reader(object):
    """Parent Reader class, to be subclassed by specific readers.
    """

    __metaclass__ = ABCMeta

    return_block = True  # By default, all readers should be
                         # cabable of returning blocks of data

    def __init__(self):
        # Common constructor for all readers

        # Set projection for coordinate transformations
        self.proj = pyproj.Proj(self.proj4)

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

    def xy2lonlat(self, x, y):
        """Calculate x,y in own projection from given lon,lat (scalars/arrays).
        """
        return self.proj(x, y, inverse=True)

    def lonlat2xy(self, lon, lat):
        """Calculate lon,lat from given x,y (scalars/arrays) in own projection.
        """
        return self.proj(lon, lat, inverse=False)


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
        if self.start_time is not None and time < self.start_time:
            raise ValueError('Requested time (%s) is before first available '
                             'time (%s) of %s' % (time, self.start_time,
                                                  self.name))
        if self.end_time is not None and time > self.end_time:
            raise ValueError('Requested time (%s) is after last available '
                             'time (%s) of %s' % (time, self.end_time,
                                                  self.name))

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

    def index_of_closest_time(self, requestedTime):
        """Return (internal) index of internal time closest to requested time.

        Assuming time step is constant; this method should be
        overloaded for readers for which this is not the case
        """
        indx = float((requestedTime - self.start_time).total_seconds()) / \
            float(self.time_step.total_seconds())
        indx = int(round(indx))
        nearestTime = self.times[indx]
        return indx, nearestTime

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
        if hasattr(self, 'depth'):
            outStr += 'Depths [m]: \n  ' + str(self.depths) + '\n'
        outStr += 'Available time range:\n'
        outStr += '  start: ' + str(self.start_time) + \
                  '   end: ' + str(self.end_time) + \
                  '   step: ' + str(self.time_step) + '\n'
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
