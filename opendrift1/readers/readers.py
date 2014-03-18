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
    pass
