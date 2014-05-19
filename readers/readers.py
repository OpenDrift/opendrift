try:
    import pyproj  # Import pyproj
except:
    try: 
        # ...alternatively use version included with Basemap
        from mpl_toolkits.basemap import pyproj
    except:
        raise ImportError('pyproj needed for coordinate transformations,'\
                          ' please install from ' \
                          'https://code.google.com/p/pyproj/')


class Readers(object):
    '''Single entry point towards driver data from several models'''

    def __init__(self):
        pass

    def set_point_reader(self, readerName, parameters):
        
        if isinstance(parameters, str):
            parameters = [parameters]
        # Import requested reader
        try:
            readerModule =__import__('reader_' + readerName)
        except ImportError:
            raise ImportError('Reader ' + readerName + ' (reader_' + 
                              readerName + '.py)' + ' not available')
            
        # Check that reader class contains the requested parameters
        reader = readerModule.Reader
        missingParameters = set(parameters) - set(reader.parameters)
        if missingParameters:
            raise ValueError('Missing parameters: ' +
                             str(list(missingParameters)))


class PointReader(object):
    '''Parent PointReader class'''

    def __init__(self):
        # Common constructor for all readers

        # Set projection for coordinate transformations
        self.proj = pyproj.Proj(self.proj4)

    def xy2lonlat(self, x, y):
       return self.proj(x, y, inverse=True)

    def lonlat2xy(self, lon, lat):
       return self.proj(lon, lat, inverse=False)

    def check_arguments(self, parameters, time, x, y, depth):
        # Check that required position and time are within
        # coverage of this reader, and that it can provide
        # the requested parameters
        for parameter in parameters:
            if parameter not in self.parameters:
                raise ValueError('Parameter not available: ' + parameter)
        if self.startTime is not None and time < self.startTime:
            raise ValueError('Outside time domain')
        if self.endTime is not None and time > self.endTime:
            raise ValueError('Outside time domain')
        if x.min() < self.xmin or x.max() > self.xmax:
            raise ValueError('x outside available domain')
        if y.min() < self.ymin or y.max() > self.ymax:
            raise ValueError('y outside available domain')

    def index_of_closest_time(self, requestedTime):
        # Temporarily assuming time steps are constant 
        indx = float((requestedTime - self.startTime).total_seconds()) / \
                float(self.timeStep.total_seconds())
        indx = int(round(indx))
        nearestTime = self.times[indx]
        return indx, nearestTime
        

    def __repr__(self):
        outStr = '===========================\n'
        outStr += 'Reader: \n'
        outStr += 'Projection: \n  ' + self.proj4 + '\n'
        outStr += '  xmin: %f   xmax: %f   step: %f\n' % (self.xmin, self.xmax, self.delta_x)
        outStr += '  ymin: %f   ymax: %f   step: %f\n' % (self.ymin, self.ymax, self.delta_y)
        corners = self.xy2lonlat([self.xmin, self.xmax, self.xmin, self.xmax],
                                 [self.ymax, self.ymax, self.ymin, self.ymin])
        outStr += '  Corners (lon, lat):\n'
        outStr += '    (%6.2f, %6.2f)  (%6.2f, %6.2f)\n' % (corners[0][0],
                                                           corners[1][0],
                                                           corners[0][2],
                                                           corners[1][2])
        outStr += '    (%6.2f, %6.2f)  (%6.2f, %6.2f)\n' % (corners[0][1],
                                                           corners[1][1],
                                                           corners[0][3],
                                                           corners[1][3]) 
        outStr += 'Available time range:\n'
        outStr += '  start: ' + str(self.startTime) + \
                    '   end: ' + str(self.endTime) + \
                    '   step: ' + str(self.timeStep) + '\n'
        outStr += 'Parameters:\n'
        for parameter in self.parameters:
            outStr += '  ' + parameter + '\n'
        outStr += '===========================\n'
        return outStr

        print self.parameters
        print self.proj4
