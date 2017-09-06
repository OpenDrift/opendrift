import numpy as np
from datetime import timedelta

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

def sensitivity_simulation(cls, lon=4.7, lat=60.0, z=0, readers=None,
                           number=1000, radius=0, seed_time=None,
                           time_step=3600, time_step_output=None,
                           duration=timedelta(hours=2),
                           filenames=None, recalculate=True):

    if recalculate is False:
        try:
            o1 = cls()
            o1.io_import_file('o0.nc')
            o2 = cls()
            o2.io_import_file('o1.nc')
            return o1, o2
        except:
            print 'Could not import'
    lon = np.atleast_1d(lon)
    if len(lon) == 1:
        lon = [lon[0], lon[0]]
    lat = np.atleast_1d(lat)
    if len(lat) == 1:
        lat = [lat[0], lat[0]]
    try:
        z[1]
    except:
        z = np.atleast_1d(z)
        if len(z) == 1:
            z = [z[0], z[0]]

    radius = np.atleast_1d(radius)

    for i in (0, 1):
        print i
        o = cls()
        o.add_readers_from_list(readers)
        if seed_time is None:
            print dir(o)
            print o.readers
            print type(o.readers)
            for r in o.readers:
                seed_time = o.readers[r].start_time
                break
        o.seed_elements(lon=lon[i], lat=lat[i], z=z[i], number=number,
                        radius=radius, time=seed_time)
        if time_step_output is None:
            time_step_output = time_step
        print o
        print duration, time_step, time_step_output
        #stop
        eargs = {'outfile': 'o%d.nc' % i}
        o.run(duration=duration, time_step=time_step,
              time_step_output=time_step_output, **eargs)
        if i == 1:
            o1 = o
        else:
            o2 = o

    return o1, o2
 
