import unittest
import importlib
import numpy as np
import time
from datetime import timedelta

from .version import __version__


# For automated access to available drift classes, e.g. for GUI
# Hardcoded for now
_available_models = \
    ['leeway.Leeway',
     'openoil.OpenOil',
     'shipdrift.ShipDrift']

def get_model_names():
    return [m.split('.')[-1] for m in _available_models]

def get_model(model_name):
    if model_name not in get_model_names():
        raise ValueError('No drift model named %s' % model_name)
    else:
        for m in _available_models:
            if m.split('.')[-1] == model_name:
                module = importlib.import_module(
                            'opendrift.models.' + m.split('.')[0])
                model = getattr(module, model_name)
                return model


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

def import_from_ladim(ladimfile, romsfile):
    """Import Ladim output file as OpenDrift simulation obejct"""

    from models.oceandrift3D import OceanDrift3D
    o = OceanDrift3D()
    from netCDF4 import Dataset, date2num, num2date
    if isinstance(romsfile, basestring):
        from opendrift.readers import reader_ROMS_native
        romsfile = reader_ROMS_native.Reader(romsfile)
    l = Dataset(ladimfile, 'r')
    pid = l.variables['pid'][:]
    particle_count = l.variables['particle_count'][:]
    end_index = np.cumsum(particle_count)
    start_index = np.concatenate(([0], end_index[:-1]))
    x = l.variables['X'][:]
    y = l.variables['Y'][:]
    lon, lat = romsfile.xy2lonlat(x, y)
    time = num2date(l.variables['time'][:],
                    l.variables['time'].units)

    history_dtype_fields = [
        (name, o.ElementType.variables[name]['dtype'])
        for name in o.ElementType.variables]
    # Add environment variables
    o.history_metadata = o.ElementType.variables.copy()
    history_dtype = np.dtype(history_dtype_fields)

    num_timesteps = len(time)
    num_elements = len(l.dimensions['particle'])
    o.history = np.ma.array(
        np.zeros([num_elements, num_timesteps]),
        dtype=history_dtype, mask=[True])

    for n in range(num_timesteps):
        start = start_index[n]
        active = pid[start:start+particle_count[n]]
        o.history['lon'][active, n] = \
            lon[start:start+particle_count[n]]
        o.history['lat'][active, n] = \
            lat[start:start+particle_count[n]]
        o.history['status'][active, n] = 0

    o.status_categories = ['active', 'missing_data']
    firstlast = np.ma.notmasked_edges(o.history['status'], axis=1)
    index_of_last = firstlast[1][1]
    o.history['status'][np.arange(len(index_of_last)),
                        index_of_last] = 1
    kwargs = {}
    for var in ['lon', 'lat', 'status']:
        kwargs[var] = o.history[var][
            np.arange(len(index_of_last)), index_of_last]
    kwargs['ID'] = range(num_elements)
    o.elements = o.ElementType(**kwargs)
    o.elements_deactivated = o.ElementType()
    o.remove_deactivated_elements()
    # Import time steps from metadata
    o.time_step = time[1] - time[0]
    o.time_step_output = o.time_step
    o.start_time = time[0]
    o.time = time[-1]
    o.steps_output = num_timesteps

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
            print('Could not import')
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
        print(i)
        o = cls()
        o.add_readers_from_list(readers)
        if seed_time is None:
            print(dir(o))
            print(o.readers)
            print(type(o.readers))
            for r in o.readers:
                seed_time = o.readers[r].start_time
                break
        o.seed_elements(lon=lon[i], lat=lat[i], z=z[i], number=number,
                        radius=radius, time=seed_time)
        if time_step_output is None:
            time_step_output = time_step
        print(o)
        print(duration, time_step, time_step_output)
        #stop
        eargs = {'outfile': 'o%d.nc' % i}
        o.run(duration=duration, time_step=time_step,
              time_step_output=time_step_output, **eargs)
        if i == 1:
            o1 = o
        else:
            o2 = o

    return o1, o2
 

# Add timer for unittest
def setUp(self):
    self._started_at = time.time()

def tearDown(self):
    elapsed = time.time() - self._started_at
    print('TIMING: ({}s) {}'.format(round(elapsed, 2), self.id()))

unittest.TestCase.setUp = setUp
unittest.TestCase.tearDown = tearDown
