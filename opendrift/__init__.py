def open(filename):
    '''Import netCDF output file as OpenDrift object of correct class'''

    import logging
    import pydoc
    from netCDF4 import Dataset
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
    o = cls()
    o.io_import_file(filename)
    logging.info('Returning ' + str(type(o)) + ' object')
    return o
