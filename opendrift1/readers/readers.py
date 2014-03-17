class Readers(object):
    '''Single entry point towards driver data from several models'''

    def __init__(object):
        pass

    def set_point_reader(self, readerName, parameters):
        
        try:
            __import__('reader_' + readerName)
        except ImportError:
            raise Error('Reader ' + readerName + ' (reader_' + readerName +
                        '.py') not available')
            
