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
        for parameter in parameters:
            if parameter not in self.parameters:
                raise ValueError('Parameter not available: ' + parameter)
        if self.startTime is not None and time < self.startTime:
            raise ValueError('Outside time domain')
        if self.endTime is not None and time < self.endTime:
            raise ValueError('Outside time domain')
        if x < self.xmin or x > self.xmax:
            raise ValueError('x outside available domain')
        if y < self.ymin or y > self.ymax:
            raise ValueError('y outside available domain')
