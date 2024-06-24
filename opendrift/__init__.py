"""
Opendrift module

.. currentmodule:: opendrift

.. doctest::

    >>> import opendrift

"""
import logging; logger = logging.getLogger(__name__)
import importlib
import numpy as np
from .version import __version__

# For automated access to available drift classes, e.g. for GUI
# Hardcoded for now
_available_models = \
    ['leeway.Leeway',
     'openoil.OpenOil',
     'larvalfish.LarvalFish',
     'plastdrift.PlastDrift',
     'shipdrift.ShipDrift',
     'openberg_old.OpenBergOld']

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


def open(filename, times=None, elements=None, load_history=True):
    '''Import netCDF output file as OpenDrift object of correct class'''

    import os
    import pydoc
    from netCDF4 import Dataset
    if not os.path.exists(filename):
        logger.info('File does not exist, trying to retrieve from URL')
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
        logger.warning(filename + ' does not contain global attributes '
                       'opendrift_module and opendrift_class, defaulting to OceanDrift')
        module_name = 'oceandrift'
        class_name = 'OceanDrift'
    n.close()

    if class_name == 'OpenOil3D':
        class_name = 'OpenOil'
        module_name = 'opendrift.models.openoil'
    if class_name == 'OceanDrift3D':
        class_name = 'OceanDrift'
        module_name = 'opendrift.models.oceandrift'
    cls = pydoc.locate(module_name + '.' + class_name)
    if cls is None:
        from opendrift.models import oceandrift
        cls = oceandrift.OceanDrift
    o = cls()
    o.io_import_file(filename, times=times, elements=elements, load_history=load_history)
    logger.info('Returning ' + str(type(o)) + ' object')
    return o

def open_xarray(filename, chunks={'trajectory': 50000, 'time': 1000}, elements=None):
    '''Import netCDF output file as OpenDrift object of correct class'''

    import os
    import pydoc
    import xarray as xr
    if not os.path.exists(filename):
        logger.info('File does not exist, trying to retrieve from URL')
        import urllib
        try:
            urllib.urlretrieve(filename, 'opendrift_tmp.nc')
            filename = 'opendrift_tmp.nc'
        except:
            raise ValueError('%s does not exist' % filename)
    n = xr.open_dataset(filename)
    try:
        module_name = n.opendrift_module
        class_name = n.opendrift_class
    except:
        raise ValueError(filename + ' does not contain '
                         'necessary global attributes '
                         'opendrift_module and opendrift_class')
    n.close()

    if class_name == 'OpenOil3D':
        class_name = 'OpenOil'
        module_name = 'opendrift.models.openoil'
    if class_name == 'OceanDrift3D':
        class_name = 'OceanDrift'
        module_name = 'opendrift.models.oceandrift'
    cls = pydoc.locate(module_name + '.' + class_name)
    if cls is None:
        from opendrift.models import oceandrift
        cls = oceandrift.OceanDrift
    o = cls()
    o.io_import_file_xarray(filename, chunks=chunks, elements=elements)

    logger.info('Returning ' + str(type(o)) + ' object')
    return o

def versions():
    import multiprocessing
    import platform
    import scipy
    import matplotlib
    import netCDF4
    import xarray
    try:
        import adios_db
        adios_version = adios_db.__version__
    except:
        adios_version = ': Not installed'
    try:
        import copernicusmarine
        copernicus_version = copernicusmarine.__version__
    except:
        copernicus_version = ': Not installed'
    import sys
    s = '\n------------------------------------------------------\n'
    s += 'Software and hardware:\n'
    s += '  OpenDrift version %s\n' % __version__
    s += '  Platform: %s, %s\n' % (platform.system(), platform.release())
    try:
        from psutil import virtual_memory
        ram = virtual_memory().total/(1024**3)
    except:
        ram = 'unknown'
    s += '  %s GB memory\n' % ram
    s += '  %s processors (%s)\n' % (multiprocessing.cpu_count(),
                                   platform.processor())
    s += '  NumPy version %s\n' % np.__version__
    s += '  SciPy version %s\n' % scipy.__version__
    s += '  Matplotlib version %s\n' % matplotlib.__version__
    s += '  NetCDF4 version %s\n' % netCDF4.__version__
    s += '  Xarray version %s\n' % xarray.__version__
    s += '  ADIOS (adios_db) version %s\n' % adios_version
    s += '  Copernicusmarine version %s\n' % copernicus_version
    s += '  Python version %s\n' % sys.version.replace('\n', '')
    s += '------------------------------------------------------\n'
    return s


def import_from_ladim(ladimfile, romsfile):
    """Import Ladim output file as OpenDrift simulation obejct"""

    from models.oceandrift import OceanDrift
    o = OceanDrift()
    from netCDF4 import Dataset, date2num, num2date
    if isinstance(romsfile, str):
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

