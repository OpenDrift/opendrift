import importlib
from abc import abstractmethod, ABCMeta
from collections import OrderedDict
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


class Readers(object):
    '''Single entry point towards driver data from several models'''

    def __init__(self):
        self.readers = {}  # Dictionary, key=name, value=reader object
        self.priority_list = OrderedDict()

    def add_reader(self, readers, variables=None):

        if isinstance(variables, str):
            variables = [variables]

        if isinstance(readers, Reader):
            readers = [readers]

        for reader in readers:
        # Check if input class is of correct type
            if not isinstance(reader, Reader):
                raise TypeError('Please provide Reader object')

            # Check that reader class contains the requested variables
            if variables is not None:
                missingVariables = set(variables) - set(reader.variables)
                if missingVariables:
                    raise ValueError('Reader %s does not provide variables: %s' \
                        % (reader.name, list(missingVariables)))

            # Finally add new reader to list
            if reader not in self.readers:
                self.readers[reader.name] = reader
                #print 'Added ' + reader.name

            # Add this reader for each of the given variables
            for variable in variables if variables else reader.variables:
                if variable in self.priority_list:
                    if reader.name not in self.priority_list[variable]:
                        self.priority_list[variable].append(reader.name)
                else:
                    self.priority_list[variable] = [reader.name]

    def get_reader_groups(self):

        # Find which groups of variables are provided by
        # the same set of readers (in the same order)
        reader_groups = []
        # Find all unique reader groups
        for variable, readers in self.priority_list.items():
            if readers not in reader_groups:
                reader_groups.append(readers)
        # Find all variables returned by the same reader group
        variable_groups = [None]*len(reader_groups)
        for variable, readers in self.priority_list.items():
            for i, readerGroup in enumerate(reader_groups):
                if readers == readerGroup:
                    if variable_groups[i]:
                        variable_groups[i].append(variable)
                    else:
                        variable_groups[i] = [variable]
        return variable_groups, reader_groups

    def list_environment_variables(self):
        variables = []
        for reader in self.readers:
            variables.extend(reader.variables)
        return variables

    def __repr__(self):

        variable_groups, reader_groups = self.get_reader_groups()
        outStr = '-------------------\n'
        outStr += 'Variables:\n'
        for i, variableGroup in enumerate(variable_groups):
            outStr += '  -----\n'
            readerGroup = reader_groups[i]
            for variable in sorted(variableGroup):
                outStr += '  ' + variable + '\n'
            for i, reader in enumerate(readerGroup):
                outStr += '     ' + str(i+1) + ') ' + reader + '\n'
        #outStr += '-------------------\n'
        return outStr

#    def get_environment(self, variables, proj4, x, y, depth, time):
#        if isinstance(variables, str):
#            variables = [variables]
#
#        missingVariables = set(variables) - set(
#            self.list_environment_variables())
#        if missingVariables:
#            raise ValueError(
#                'Please add readers for the following variables: ' +
#                str(list(missingVariables)))
#
#        # Use first reader - TEMPORARY SOLUTION
#        stop
#        env = self.readers[0].get_variables(variables, time, x, y, depth)
#        #print env
#        #stop

class Reader(object):
    __metaclass__ = ABCMeta
    '''Parent Reader class'''

    def __init__(self):
        # Common constructor for all readers

        # Set projection for coordinate transformations
        self.proj = pyproj.Proj(self.proj4)

    @abstractmethod
    def get_variables(self, requestedVariables, time=None,
                      x=None, y=None, depth=None):
        pass

    def xy2lonlat(self, x, y):
        return self.proj(x, y, inverse=True)

    def lonlat2xy(self, lon, lat):
        return self.proj(lon, lat, inverse=False)

    def check_arguments(self, variables, time, x, y, depth):
        # Check that required position and time are within
        # coverage of this reader, and that it can provide
        # the requested variable(s)

        # Check time
        if time is None:
            time = self.startTime  # Get data from first timestep, if not given
            indxTime = 0

        # Convert variables to list and x,y to ndarrays
        if isinstance(variables, str):
            variables = [variables]
        x = np.asarray(x)
        y = np.asarray(y)
        depth = np.asarray(depth)

        for variable in variables:
            if variable not in self.variables:
                raise ValueError('Variable not available: ' + variable)
        if self.startTime is not None and time < self.startTime:
            raise ValueError('Requested time (%s) is before first available '
                             'time (%s)' % (time, self.startTime))
        if self.endTime is not None and time > self.endTime:
            raise ValueError('Requested time (%s) is after last available '
                             'time (%s)' % (time, self.endTime))
        #if x.min() < self.xmin or x.max() > self.xmax:
        #    raise ValueError('Requested x <%s, %s> is outside domain '
        #            '<%s, %s>' % (x.min(), x.max(), self.xmin, self.xmax))
        #if y.min() < self.ymin or y.max() > self.ymax:
        #    raise ValueError('Requested y <%s, %s> is outside domain '
        #            '<%s, %s>' % (y.min(), y.max(), self.ymin, self.ymax))
        return variables, time, x, y, depth

    def index_of_closest_time(self, requestedTime):
        # Assuming time steps are constant; this method should be
        # overloaded for readers for which this is not the case
        indx = float((requestedTime - self.startTime).total_seconds()) / \
            float(self.timeStep.total_seconds())
        indx = int(round(indx))
        nearestTime = self.times[indx]
        return indx, nearestTime

    def indices_min_max_depth(self, depths):
        # Return min and max indices of depth dimension,
        # covering the requested depths
        minIndex = (self.depths <= depths.min()).argmin() - 1
        maxIndex = (self.depths >= depths.max()).argmax()
        return minIndex, maxIndex

    def __repr__(self):
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
        outStr += '  start: ' + str(self.startTime) + \
                  '   end: ' + str(self.endTime) + \
                  '   step: ' + str(self.timeStep) + '\n'
        outStr += 'Variables:\n'
        for variable in self.variables:
            outStr += '  ' + variable + '\n'
        outStr += '===========================\n'
        return outStr
