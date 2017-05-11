def open(filename):
    '''Import netCDF output file as OpenDrift object of correct class'''

    import os
    import logging
    import pydoc
    from netCDF4 import Dataset
    if not os.path.exists(filename):
        logging.info('File does not exist, trying to retrieve from URL')
        import urllib
        try:
            urllib.urlretrieve(filename, 'opendrift_tmp.nc')
            filename = 'opendrift_tmp.nc'
        except:
            raise ValueError('%s does not exist' % filename)
    n = Dataset(filename)
    try:
        module_name = n.opendrift_module
        class_name = n.opendrift_class
    except:
        raise ValueError(filename + ' does not contain '
                         'necessary global attributes '
                         'opendrift_module and opendrift_class')
    n.close()

    cls = pydoc.locate(module_name + '.' + class_name)
    if cls is None:
        from models import oceandrift3D
        cls = oceandrift3D.OceanDrift3D
    o = cls()
    o.io_import_file(filename)
    logging.info('Returning ' + str(type(o)) + ' object')
    return o
