# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway
# Copyright 2020, Gaute Hope, MET Norway

import sys
import os
import glob
import types
import traceback
import inspect
import logging; logging.captureWarnings(True)
import warnings
from datetime import datetime, timedelta
from collections import OrderedDict
from abc import ABCMeta, abstractmethod, abstractproperty
import netCDF4
import nc_time_axis
import xarray as xr
from netCDF4 import Dataset, date2num

import numpy as np
import scipy
import pyproj

import opendrift
from opendrift.readers.basereader import BaseReader, vector_pairs_xy
from opendrift.readers import reader_from_url
from opendrift.models.physics_methods import PhysicsMethods
from opendrift.analysis.plot import Plotter, Animator

from .seed import Seeder

class OpenDriftSimulation(PhysicsMethods, Seeder, Animator, Plotter):
    """Generic trajectory model class, to be extended (subclassed).

    This as an Abstract Base Class, meaning that only subclasses can
    be initiated and used.
    Any specific subclass ('model') must contain its own (or shared)
    specific type of particles (ElementType), whose properties are
    updated at each time_step using method update() on basis of model
    physics/chemistry/biology and 'required_variables' (environment)
    which are provided by one or more Reader objects.

    Attributes:
        ElementType: the type (class) of particles to be used by this model

        elements: object of the class ElementType, storing the specific
        particle properties (ndarrays and scalars) of all active particles
        as named attributes. Elements are added by seeding-functions
        (presently only one implemented: seed_elements).

        elements_deactivated: ElementType object containing particles which
            have been deactivated (and removed from 'elements')

        elements_scheduled: ElementType object containing particles which
            have been scheduled, but not yet activated
        required_variables: list of strings of CF standard_names which is
            needed by this model (update function) to update properties of
            particles ('elements') at each time_step. This core class has
            no required_elements, this is implemented by subclasses/modules.
        environment: recarray storing environment variables (wind, waves,
            current etc) as named attributes. Attribute names follow
            standard_name from CF-convention, allowing any OpenDriftSimulation
            module/subclass using environment data from any readers which
            can provide the requested variables. Used in method 'update'
            to update properties of elements every time_step.
        proj4: string defining the common spatial reference system (SRS) onto
            which data from all readers are interpolated
        proj: Proj object initialised from proj4 string; used for
            coordinate tranformations
        time_step: timedelta object, time interval at which element properties
            are updated (including advection).
        time_step_output: timedelta object, time interval at which element
            properties are stored in memory and eventually written to file
        readers: Dictionary where values are Reader objects, and names are
            unique reference keywords used to access a given reader (typically
            filename or URL)
        priority_list: OrderedDict where names are variable names,
            and values are lists of names (kewywords) of the reader, in the
            order of priority (user defined) of which the readers shall be
            called to retrieve environmental data.
    """

    __metaclass__ = ABCMeta

    status_categories = ['active']  # Particles are active by default

    # Default plotting colors of trajectory endpoints
    status_colors_default = {'initial': 'green',
                             'active': 'blue',
                             'missing_data': 'gray'}

    CONFIG_LEVEL_ESSENTIAL=1
    CONFIG_LEVEL_BASIC=2
    CONFIG_LEVEL_ADVANCED=3

    max_speed = 1  # Assumed max average speed of any element
    required_profiles_z_range = None  # [min_depth, max_depth]
    plot_comparison_colors = ['k', 'r', 'g', 'b', 'm', 'c', 'y']

    logger = logging.getLogger('opendrift')

    def __init__(self, proj4=None, seed=0, iomodule='netcdf',
                 loglevel=logging.DEBUG, logtime='%H:%M:%S', logfile=None):
        """Initialise OpenDriftSimulation

        Args:
            proj4: proj4 string defining spatial reference system.
                If not specified, SRS is taken from the first added reader.
            seed: integer or None. A given integer will yield identical
                random numbers drawn each simulation. Random numbers are
                e.g. used to distribute particles spatially when seeding,
                and may be used by modules (subclasses) for e.g. diffusion.
                Specifying a fixed value (default: 0) is useful for sensitivity
                tests. With seed = None, different random numbers will be drawn
                for subsequent runs, even with identical configuration/input.
            iomodule: name of module used to export data
                default: netcdf, see folder 'io' for more alternatives.
                'iomodule' is module/filename without preceeding 'io_'
            loglevel: set to 0 (default) to retrieve all debug information.
                Provide a higher value (e.g. 20) to receive less output.
                Use the string 'custom' to configure logging from outside.
            logtime: if True, a time stamp is given for each logging line.
                     logtime can also be given as a python time specifier
                     (e.g. '%H:%M:%S')
        """

        self.show_continuous_performance = False

        # Dict to store readers
        self.readers = OrderedDict()  # Dictionary, key=name, value=reader object
        self.priority_list = OrderedDict()
        self.use_block = True  # Set to False if interpolation left to reader

        # Make copies of dictionaries so that they are private to each instance
        self.status_categories = ['active']  # Particles are active by default
        self.status_colors_default = self.status_colors_default.copy()

        if hasattr(self, 'status_colors'):
            # Append model specific colors to (and override) default colors
            self.status_colors_default.update(self.status_colors)
            self.status_colors = self.status_colors_default
        else:
            self.status_colors = self.status_colors_default

        # Set projection, if given
        self.set_projection(proj4)

        # Using a fixed seed will generate the same random numbers
        # each run, useful for sensitivity tests
        # Use seed = None to get different random numbers each time
        np.random.seed(seed)

        self.steps_calculation = 0  # Increase for each simulation step
        self.steps_output = 0
        self.elements_deactivated = self.ElementType()  # Empty array
        self.elements = self.ElementType()  # Empty array

        if loglevel != 'custom':
            if logfile is not None:
                handler = logging.FileHandler(logfile, mode='w')
            else:
                handler = logging.StreamHandler()
            format = '%(levelname)s: %(message)s'
            datefmt = None
            if logtime is not False:
                format = '%(asctime)s ' + format
                if logtime is not True:
                    datefmt = logtime
            formatter = logging.Formatter(format, datefmt=datefmt)
            handler.setFormatter(formatter)
            if loglevel < 10:  # 0 is NOTSET, giving no output
                loglevel=10
            self.logger.setLevel(loglevel)
            self.logger.handlers = []
            self.logger.addHandler(handler)
            self.logger.propagate = False

        # Prepare outfile
        try:
            io_module = __import__('opendrift.export.io_' + iomodule,
                                   fromlist=['init', 'write_buffer',
                                             'close', 'import_file'])
        except ImportError:
            self.logger.info('Could not import iomodule ' + iomodule)
        self.io_init = types.MethodType(io_module.init, self)
        self.io_write_buffer = types.MethodType(io_module.write_buffer, self)
        self.io_close = types.MethodType(io_module.close, self)
        self.io_import_file = types.MethodType(io_module.import_file, self)
        self.io_import_file_xarray = types.MethodType(io_module.import_file_xarray, self)

        # Set configuration options
        self._add_config({
            # type, default, min, max, enum, important, value, units, description
            'general:use_auto_landmask': {'type': 'bool', 'default': True,
                'description': 'A built-in GSHHG global landmask is used if True, '
                 'otherwise landmask is taken from reader or fallback value.',
                 'level': self.CONFIG_LEVEL_ADVANCED},
            'general:coastline_action': {'type': 'enum', 'enum': ['none', 'stranding', 'previous'],
                'default': 'stranding', 'level': self.CONFIG_LEVEL_BASIC,
                 'description': 'None means that objects may also move over land. '
                    'stranding means that objects are deactivated if they hit land. '
                    'previous means that objects will move back to the previous location '
                    'if they hit land'},
            'general:time_step_minutes': {'type': 'float', 'min': .01, 'max': 1440, 'default': 60,
                'units': 'minutes', 'level': self.CONFIG_LEVEL_BASIC, 'description':
                'Calculation time step used for the simulation. The output time step may '
                'be equal or larger than this.'},
            'general:time_step_output_minutes': {'type': 'float', 'min': 1, 'max': 1440, 'default': None,
                'units': 'minutes', 'level': self.CONFIG_LEVEL_BASIC, 'description':
                'Output time step, i.e. the interval at which output is saved. This must be larger than '
                'the calculation time step, and be an integer multiple of this.'},
            'seed:ocean_only': {'type': 'bool', 'default': True,
                'description': 'If True, elements seeded on land will be moved to the closest '
                    'position in ocean.', 'level': self.CONFIG_LEVEL_ADVANCED},
            'seed:number': {'type': 'int', 'default': 1,
                'min': 1, 'max': 100000000, 'units': 1,
                'description': 'The number of elements for the simulation.',
                'level': self.CONFIG_LEVEL_BASIC},
            'drift:max_age_seconds': {'type': 'float', 'default': None,
                'min': 0, 'max': np.inf, 'units': 'seconds',
                'description': 'Elements will be deactivated when this age is reached..',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:scheme': {'type': 'enum', 'enum': ['euler', 'runge-kutta', 'runge-kutta4'], 'default': 'euler',
                'level': self.CONFIG_LEVEL_ADVANCED, 'description':
                'Numerical advection scheme for ocean current advection.'},
            'drift:current_uncertainty': {'type': 'float', 'default': 0,
                'min': 0, 'max': 5, 'units': 'm/s',
                'description': 'Add gaussian perturbation with this standard deviation to current components at each time step.',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:current_uncertainty_uniform': {'type': 'float', 'default': 0,
                'min': 0, 'max': 5, 'units': 'm/s',
                'description': 'Add gaussian perturbation with this standard deviation to current components at each time step.',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:horizontal_diffusivity': {'type': 'float', 'default': 0,
                'min': 0, 'max': 100, 'units': 'm2/s',
                'description': 'Add horizontal diffusivity (random walk)',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:wind_uncertainty': {'type': 'float', 'default': 0,
                'min': 0, 'max': 5, 'units': 'm/s',
                'description': 'Add gaussian perturbation with this standard deviation to wind components at each time step.',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:relative_wind': {'type': 'bool', 'default': False,
                'description': 'If True, wind drift is calculated for absolute wind (wind vector minus ocean surface current vector).',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:deactivate_north_of': {'type': 'float', 'default': None,
                'min': -90, 'max': 90, 'units': 'degrees',
                'description': 'Elements are deactivated if the move further north than this limit',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:deactivate_south_of': {'type': 'float', 'default': None,
                'min': -90, 'max': 90, 'units': 'degrees',
                'description': 'Elements are deactivated if the move further south than this limit',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:deactivate_east_of': {'type': 'float', 'default': None,
                'min': -360, 'max': 360, 'units': 'degrees',
                'description': 'Elements are deactivated if the move further east than this limit',
                'level': self.CONFIG_LEVEL_ADVANCED},
            'drift:deactivate_west_of': {'type': 'float', 'default': None,
                'min': -360, 'max': 360, 'units': 'degrees',
                'description': 'Elements are deactivated if the move further west than this limit',
                'level': self.CONFIG_LEVEL_ADVANCED},
        })

        # Add default element properties to config
        c = {}
        for p in self.ElementType.variables:
            v = self.ElementType.variables[p]
            if 'seed' in v and v['seed'] is False:
                continue  # Properties which may not be provided by user
            minval = v['min'] if 'min' in v else None
            maxval = v['max'] if 'max' in v else None
            units = v['units'] if 'units' in v else None
            c['seed:%s' % p] = {
                'type': v['type'] if 'type' in v else 'float',
                'min': v['min'] if 'min' in v else None,
                'max': v['max'] if 'max' in v else None,
                'units': v['units'] if 'units' in v else None,
                'default': v['default'] if 'default' in v else None,
                'description': v['description'] if 'description' in v else 'Seeding value of %s' % p,
                'level': v['level'] if 'level' in v else self.CONFIG_LEVEL_ADVANCED}
        self._add_config(c)

        # Add constant and fallback environment variables to config
        c = {}
        for v in self.required_variables:
            c['environment:constant:%s' % v] = {'type': 'float',
                'min': None, 'max': None, 'units': None, 'default': None,
                'level': OpenDriftSimulation.CONFIG_LEVEL_BASIC,
                'description': 'Use constant value for %s' %v}
            c['environment:fallback:%s' % v] = {'type': 'float',
                'min': None, 'max': None, 'units': None, 'default':
                 self.required_variables[v]['fallback'] if
                    'fallback' in self.required_variables[v] else None,
                'level': OpenDriftSimulation.CONFIG_LEVEL_BASIC,
                'description': 'Fallback value for %s if not available from any reader' %v}
        self._add_config(c)

        # Find variables which require profiles
        self.required_profiles = [var for var in self.required_variables
            if 'profiles' in self.required_variables[var] and
            self.required_variables[var]['profiles'] is True]

        # Find variables which are desired, but not required
        self.desired_variables = [var for var in self.required_variables
            if 'important' in self.required_variables[var] and
            self.required_variables[var]['important'] is False]

        self.timer_start('total time')
        self.timer_start('configuration')

        self.add_metadata('opendrift_version', opendrift.__version__)
        self.logger.info('OpenDriftSimulation initialised (version %s)' %
                     opendrift.__version__)

    def list_config(self, prefix=''):
        """List all possible configuration settings with values"""
        str = '\n=============================================\n'
        for key in self._config:
            if key.startswith(prefix):
                str += '%s [%s]\n' % (key, self.get_config(key))
        str += '=============================================\n'
        self.logger.info(str)

    def list_configspec(self, prefix=''):
        """Readable formatting of config specification"""
        for c, i in self._config.items():
            if c.startswith(prefix):
                val = i['value'] if 'value' in i else None
                if i['type'] == 'bool':
                    rang=''
                elif i['type'] in ['float', 'int']:
                    rang = 'min: %s, max: %s [%s]' % (i['min'], i['max'], i['units'])
                elif i['type'] == 'enum':
                    rang = i['enum']
                print('%-35s [%s] %-5s %s %s...' % (c, val, i['type'], rang, i['description'][0:20]))

    def get_configspec(self, prefix='', level=[1,2,3]):
        if not isinstance(level, list):
            level = [level]
        configspec = {k:v for (k,v) in self._config.items() if k.startswith(prefix)
                        and self._config[k]['level'] in level}
        return configspec

    def _add_config(self, config, overwrite=True):
        """Add configuration settings

        config is a dictionary where keys are configuration keywords,
        and values are dictionaries with the following contents:

        type (string): 'float', 'int', 'bool' or 'enum'

        min, max (float/int/None): (only when type is 'float' or 'int')
            The minimum and maximum allowed values for this setting.
            May also be None if there are no upper/lowe limits.

        units (string): (only when type is 'float' or 'int')
            The units of this config setting.

        enum (list): (only when type is 'enum')
            A list of possible values for this setting.

        default (number/bool/string/None):
            The default value for this setting.

        value (number/bool/string/None): The actual value for this setting.
            This is updated with self.set_config(key, value) and retrieved
            with self.get_config(key)

        description (string):
            A description of this config setting, for users/documentation/GUIs.

        level (int): A parameter to determine the level of exposure in GUIs
            1 self.CONFIG_LEVEL_ESSENTIAL: important setting which user has to consider
            2 self.CONFIG_LEVEL_BASIC: setting which many users may consider
            3 self.CONFIG_LEVEL_ADVANCED: setting relevant only to advanced users

        """

        caller = inspect.stack()[1]
        caller = os.path.splitext(os.path.basename(caller.filename))[0]
        self.logger.debug('Adding %i config items from %s' % (len(config), caller))
        if not hasattr(self, '_config'):
            self._config = {}
        remove = []
        for c, i in config.items():  # Check that provided config is conistent
            if c in self._config:
                if overwrite is False:
                    self.logger.debug('  Config item %s is already specified, not overwriting' % c)
                    remove.append(c)
                else:
                    self.logger.debug('  Overwriting config item %s' % c)
            for p in ['type', 'description', 'level']:
                if p not in i:
                    raise ValueError('"%s" must be specified for config item %s' % (p, c))
            if i['level'] != self.CONFIG_LEVEL_ESSENTIAL and 'default' not in i:#or i['default'] is None:
                raise ValueError('A default value must be provided for config item %s' % c)
            if i['type'] == 'enum':
                if 'enum' not in i or not isinstance(i['enum'], list):
                    raise ValueError('"enum" of type list must be provided for config item %s' % (c))
            elif i['type'] in ['float', 'int']:
                for p in ['min', 'max', 'units']:
                    if p not in i:
                        raise ValueError('"%s" not provided for config item %s' % (p, c))
            elif i['type'] == 'bool':
                pass  # no check for bool
            else:
                raise ValueError('Config type "%s" (%s) is not defined. Valid options are: '
                                 'float, int, enum, bool' % (i['type'], c))
            if 'default' in i:
                i['value'] = i['default']
        for r in remove:
            del config[r]
        self._config.update(config)

    def set_config(self, key, value):
        if not key in self._config:
            raise ValueError('No config setting named %s' % key)
        i = self._config[key]
        if i['type'] == 'bool':
            if value not in [True, False]:
                raise ValueError('Config value %s must be True or False' % key)
        elif i['type'] in ['float', 'int']:
            if (i['min'] is not None and value < i['min']) or (i['max'] is not None and value > i['max']):
                raise ValueError('Config value %s must be between %s and %s' % (key, i['min'], i['max']))
            if i['type'] == 'float' and value is not None:
                value = np.float(value)
            elif i['type'] == 'int' and value is not None:
                value = np.int(value)
        elif i['type'] == 'enum':
            if value not in i['enum']:
                if len(i['enum']) > 5:
                    import difflib
                    matches = difflib.get_close_matches(value, i['enum'], n=20, cutoff=.3)
                    containing = [e for e in i['enum'] if value in e]
                    matches = list(set(matches) | set(containing))
                    if len(matches) > 0:
                        matches.sort()
                        suggestion = '\nDid you mean any of these?\n%s' % str(matches)
                else:
                    suggestion = ''
                raise ValueError('Wrong configuration, possible values are:\n\t%s\n%s' %
                                 (i['enum'], suggestion))

        self._config[key]['value'] = value

    def _set_config_default(self, key, value):
        """Update both default and actual value of a config setting"""
        self.set_config(key, value)
        self._config[key]['default'] = self.get_config(key)

    def get_config(self, key):
        if not key in self._config:
            raise ValueError('No config setting named %s' % key)
        return(self._config[key]['value'])

    def add_metadata(self, key, value):
        """Add item to metadata dictionary, for export as netCDF global attributes"""
        if not hasattr(self, 'metadata_dict'):
            from collections import OrderedDict
            self.metadata_dict = OrderedDict()
        self.metadata_dict[key] = value

    def prepare_run(self):
        pass  # to be overloaded when needed

    def store_present_positions(self, IDs=None, lons=None, lats=None):
        """Store present element positions, in case they shall be moved back"""
        if self.get_config('general:coastline_action') == 'previous':
            if not hasattr(self, 'previous_lon'):
                self.previous_lon = np.ma.masked_all(self.num_elements_total())
                self.previous_lat = np.ma.masked_all(self.num_elements_total())
            if IDs is None:
                IDs = self.elements.ID
                lons = self.elements.lon
                lats = self.elements.lat
                self.newly_seeded_IDs = None
            else:
                # to check if seeded on land
                if len (IDs) > 0:
                    self.newly_seeded_IDs = np.copy(IDs)
                else:
                    self.newly_seeded_IDs = None
            self.previous_lon[IDs-1] = np.copy(lons)
            self.previous_lat[IDs-1] = np.copy(lats)

    def store_previous_variables(self):
        """Store some environment variables, for access at next time step"""

        if not hasattr(self, 'store_previous'):
            return
        if not hasattr(self, 'variables_previous'):
            # Create ndarray to store previous variables
            dtype = [(var, np.float32) for var in self.store_previous]
            self.variables_previous = np.array(
                np.full(self.num_elements_total(), np.nan), dtype=dtype)

        # Copying variables_previous to environment_previous
        self.environment_previous = self.variables_previous[self.elements.ID-1]

        # Use new values for new elements which have no previous value
        for var in self.store_previous:
            undefined = np.isnan(self.environment_previous[var])
            self.environment_previous[var][undefined] = getattr(self.environment, var)[undefined]

        self.environment_previous = self.environment_previous.view(np.recarray)

        for var in self.store_previous:
            self.variables_previous[var][self.elements.ID-1] = getattr(self.environment, var)

    def interact_with_coastline(self, final=False):
        if self.num_elements_active() == 0:
            return
        """Coastline interaction according to configuration setting"""
        i = self.get_config('general:coastline_action')
        if not hasattr(self, 'environment') or not hasattr(self.environment, 'land_binary_mask'):
            return
        if i == 'none':  # Do nothing
            return
        if final is True:  # Get land_binary_mask for final location
            en, en_prof, missing = \
                self.get_environment(['land_binary_mask'],
                                     self.time,
                                     self.elements.lon,
                                     self.elements.lat,
                                     self.elements.z,
                                     None)
            self.environment.land_binary_mask = en.land_binary_mask

        if i == 'stranding':  # Deactivate elements on land
            self.deactivate_elements(
                self.environment.land_binary_mask == 1, reason='stranded')
        elif i == 'previous':  # Go back to previous position (in water)
            if self.newly_seeded_IDs is not None:
                self.deactivate_elements(
                    (self.environment.land_binary_mask == 1) &
                    (self.elements.age_seconds == self.time_step.total_seconds()),
                    reason='seeded_on_land')
            on_land = np.where(self.environment.land_binary_mask == 1)[0]
            if len(on_land) == 0:
                self.logger.debug('No elements hit coastline.')
            else:
                self.logger.debug('%s elements hit coastline, '
                              'moving back to water' % len(on_land))
                on_land_ID = self.elements.ID[on_land]
                self.elements.lon[on_land] = \
                    np.copy(self.previous_lon[on_land_ID - 1])
                self.elements.lat[on_land] = \
                    np.copy(self.previous_lat[on_land_ID - 1])
                self.environment.land_binary_mask[on_land] = 0

    @abstractmethod
    def update(self):
        """Any trajectory model implementation must define an update method.
        This method must/can use environment data (self.environment) to
        update properties (including position) of its particles (self.elements)
        """

    @abstractproperty
    def ElementType(self):
        """Any trajectory model implementation must define an ElementType."""

    @abstractproperty
    def required_variables(self):
        """Any trajectory model implementation must list needed variables."""

    def test_data_folder(self):
        import opendrift
        return os.path.abspath(
            os.path.join(os.path.dirname(opendrift.__file__),
                         '..', 'tests', 'test_data')) + os.path.sep

    def set_projection(self, proj4):
        """Set the projection onto which data from readers is reprojected."""
        self.proj4 = proj4
        if proj4 is not None:
            self.proj = pyproj.Proj(self.proj4 + ' +ellps=WGS84')
            self.logger.debug('Calculation SRS set to: ' + self.proj.srs)
        else:
            self.proj = None
            self.logger.debug('Calculation SRS set to: ' + str(self.proj))

    def lonlat2xy(self, lon, lat):
        """Calculate x,y in own projection from given lon,lat (scalars/arrays).
        """
        if 'ob_tran' in self.proj4:
            x, y = self.proj(lon, lat, inverse=False)
            return np.degrees(x), np.degrees(y)
        elif self.proj.crs.is_geographic:
            return lon, lat
        else:
            x, y = self.proj(lon, lat, inverse=False)
            return x, y

    def xy2lonlat(self, x, y):
        """Calculate lon,lat from given x,y (scalars/arrays) in own projection.
        """
        if self.proj.crs.is_geographic:
            if 'ob_tran' in self.proj4:
                self.logger.info('NB: Converting deg to rad due to ob_tran srs')
                x = np.radians(np.array(x))
                y = np.radians(np.array(y))
                return self.proj(x, y, inverse=True)
            else:
                return x, y
        else:
            return self.proj(x, y, inverse=True)

    def timer_start(self, category):
        if not hasattr(self, 'timers'):
            self.timers = OrderedDict()
        if not hasattr(self, 'timing'):
            self.timing = OrderedDict()
        if category not in self.timing:
            self.timing[category] = timedelta(0)
        self.timers[category] = datetime.now()

    def timer_end(self, category):
        if self.timers[category] is not None:
            self.timing[category] += datetime.now() - self.timers[category]
        self.timers[category] = None

    def format_timedelta(self, timedelta):
        '''Format timedelta nicely for display'''
        timestr = str(timedelta)[0:str(timedelta).find('.') + 2]
        for i, c in enumerate(timestr):
            if c in '123456789.':
                timestr = timestr[i:]  # Strip leading 0 and :
                if c == '.':
                    timestr = '0' + timestr
                break
        return timestr

    def performance(self):
        '''Report the time spent on various tasks'''

        outStr = '--------------------\n'
        outStr += 'Reader performance:\n'
        for r in self.readers:
            reader = self.readers[r]
            if reader.is_lazy:
                continue
            outStr += '--------------------\n'
            outStr += r + '\n'
            outStr += reader.performance()

        outStr += '--------------------\n'
        outStr += 'Performance:\n'
        for category, time in self.timing.items():
            timestr = str(time)[0:str(time).find('.') + 2]
            for i, c in enumerate(timestr):
                if c in '123456789.':
                    timestr = timestr[i:]  # Strip leading 0 and :
                    if c == '.':
                        timestr = '0' + timestr
                    break
            parts = category.split(':')
            indent = '  '*(len(parts) - 1)
            category = parts[-1]
            category = category.replace('<colon>', ':')
            outStr += '%s%7s %s\n' % (indent, timestr, category)

        outStr += '--------------------\n'
        return outStr

    def add_reader(self, readers, variables=None, first=False):
        """Add one or more readers providing variables used by this model.

        Method may be called subsequently to add more readers
        for other variables.

        Args:
            readers: one or more (list) Reader objects.
            variables: optional, list of strings of standard_name of
                variables to be provided by this/these reader(s).
            first: Set to True if this reader should be set as first option
        """

        # Convert any strings to lists, for looping
        if isinstance(variables, str):
            variables = [variables]
        if isinstance(readers, BaseReader):
            readers = [readers]

        for reader in readers:
            # Check if input class is of correct type
            if not isinstance(reader, BaseReader) and \
                    not hasattr(reader, '_lazyname'):
                raise TypeError('Please provide Reader object')

            # Check that reader class contains the requested variables
            if variables is not None:
                missingVariables = set(variables) - set(reader.variables)
                if missingVariables:
                    raise ValueError('Reader %s does not provide variables: %s'
                                     % (reader.name, list(missingVariables)))

            # Finally add new reader to list
            if reader.name in self.readers:
                # Reader names must be unique, adding integer
                for n in range(99999):
                    tmp_name = reader.name + '_%d' % n
                    if tmp_name not in self.readers:
                        reader.name = tmp_name
                        break

            # Horizontal buffer of reader must be large enough to cover
            # the distance possibly covered by elements within a time step
            if not reader.is_lazy:
                reader.set_buffer_size(max_speed=self.max_speed)

            self.readers[reader.name] = reader
            if self.proj is None and not reader.is_lazy:
                if reader.proj4 is not None and reader.proj4 != 'None':
                    self.set_projection(reader.proj4)
                    self.logger.debug('Using srs for common grid: %s' %
                                  self.proj4)
                else:
                    self.logger.debug('%s is unprojected, cannot use '
                                  'for common grid' % reader.name)
            self.logger.debug('Added reader ' + reader.name)

            # Add this reader for each of the given variables
            if reader.is_lazy is False:
                for variable in variables if variables else reader.variables:
                    if variable in list(self.priority_list):
                        if reader.name not in self.priority_list[variable]:
                            if first is True:
                                self.priority_list[variable].insert(0, reader.name)
                            else:
                                self.priority_list[variable].append(reader.name)
                    else:
                        self.priority_list[variable] = [reader.name]

        # Remove/hide variables not needed by the current trajectory model
        for variable in list(self.priority_list):
            if variable not in self.required_variables:
                del self.priority_list[variable]

        # Set projection to latlong if not taken from any of the readers
        if self.proj is None:
            self.logger.info('Setting SRS to latlong, since not defined before.')
            self.set_projection('+proj=latlong')

    def add_readers_from_list(self, urls, timeout=10, lazy=True):
        '''Make readers from a list of URLs or paths to netCDF datasets'''

        if isinstance(urls, str):
            urls = [urls]
        if lazy is True:
            from opendrift.readers.reader_lazy import Reader
            readers = [Reader(u) for u in urls]
            self.add_reader(readers)
            return

        readers = [reader_from_url(u, timeout) for u in urls]
        self.add_reader([r for r in readers if r is not None])

    def add_readers_from_file(self, filename, timeout=10, lazy=True):
        fp = open(filename, 'r')
        sources = fp.readlines()
        sources = [line.strip() for line in sources if line[0] != '#']
        self.add_readers_from_list(sources, timeout, lazy=lazy)

    def list_environment_variables(self):
        """Return list of all variables provided by the added readers."""
        variables = []
        for reader in self.readers:
            variables.extend(self.readers[reader].variables)
        return variables

    def missing_variables(self):
        """Return list of all variables for which no reader has been added."""
        return [var for var in self.required_variables
                if var not in self.priority_list]

    def get_reader_groups(self, variables=None):
        """Find which groups of variables are provided by the same readers.

        This function loops through 'priority_list' (see above) and groups
        all variables returned by the same readers in the same order. This
        allows asking readers for several variables simultaneously,
        improving performance. Used by method 'get_environment'.

        Returns:
            variable_groups: list of lists of (environment) variables.
            reader_groups: list of list of reader names, corresponding to
            each of the variable_groups.
        """
        if variables is None:
            variables = list(self.required_variables)
        reader_groups = []
        # Find all unique reader groups
        for variable, readers in self.priority_list.items():
            if (variable in variables) and (readers not in reader_groups):
                reader_groups.append(readers)
        # Find all variables returned by the same reader group
        variable_groups = [None]*len(reader_groups)
        for variable, readers in self.priority_list.items():
            if variable not in variables:
                continue
            for i, readerGroup in enumerate(reader_groups):
                if readers == readerGroup:
                    if variable_groups[i]:
                        variable_groups[i].append(variable)
                    else:
                        variable_groups[i] = [variable]

        missing_variables = list(set(variables) -
                                 set(self.priority_list.keys()))

        return variable_groups, reader_groups, missing_variables

    def _lazy_readers(self):
        return [r for r in self.readers
                if self.readers[r].is_lazy is True]

    def _unlazy_readers(self):
        return [r for r in self.readers
                if self.readers[r].is_lazy is False]

    def _initialise_next_lazy_reader(self):
        '''Returns reader if successful and None if no more readers'''

        lazy_readers = self._lazy_readers()

        if len(lazy_readers) == 0:
            return None

        lazyname = lazy_readers[0]
        reader = self.readers[lazyname]

        try:
            reader.initialise()
        except Exception as e:
            self.logger.debug(e)
            self.logger.warning('Reader could not be initialised, and is'
                            ' discarded: ' + lazyname)
            self.discard_reader(reader)
            return self._initialise_next_lazy_reader()  # Call self

        reader.set_buffer_size(max_speed=self.max_speed)
        # Update reader lazy name with actual name
        self.readers[reader.name] = \
            self.readers.pop(lazyname)
        for var in reader.variables:
            if var in list(self.priority_list):
                self.priority_list[var].append(reader.name)
            else:
                self.priority_list[var] = [reader.name]
        # Remove variables not needed
        for variable in list(self.priority_list):
            if variable not in self.required_variables:
                del self.priority_list[variable]

        return reader

    def earliest_time(self):
        return min(self.start_time, self.expected_end_time)

    def latest_time(self):
        return min(self.start_time, self.expected_end_time)

    def discard_reader_if_not_relevant(self, reader):
        if hasattr(self, 'expected_endtime') and reader.start_time is not None and (
                (reader.start_time > self.latest_time() or
                 reader.end_time < self.earliest_time())):
            self.logger.debug('Reader does not cover simulation period')
            self.discard_reader(reader)
            return True
        if len(set(self.required_variables) &
               set(reader.variables)) == 0:
            self.logger.debug('Reader does not contain any relevant variables')
            self.discard_reader(reader)
            return True
        return False

    def discard_reader(self, reader):
        readername = reader.name
        self.logger.debug('Discarding reader: ' + readername)
        del self.readers[readername]
        if not hasattr(self, 'discarded_readers'):
            self.discarded_readers = [readername]
        else:
            self.discarded_readers.append(readername)

        # Remove from priority list
        for var in self.priority_list:
            self.priority_list[var] = [r for r in
                self.priority_list[var] if r != readername]
            if len(self.priority_list[var]) == 0:
                del self.priority_list[var]

    def discard_irrelevant_readers(self):
        for readername in self.readers:
            reader = self.readers[readername]
            if reader.is_lazy:
                continue
            if self.discard_reader_if_not_relevant(reader):
                self.logger.debug('DISCARDED: ' + readername)

    def get_environment(self, variables, time, lon, lat, z, profiles):
        '''Retrieve environmental variables at requested positions.

        Updates:
            Buffer (raw data blocks) for each reader stored for performance:
                [readers].var_block_before (last before requested time)
                [readers].var_block_after (first after requested time)
                - lists of one ReaderBlock per variable group: time, x, y, [vars]

        Returns:
            environment: recarray with variables as named attributes,
                         interpolated to requested positions/time.

        '''
        self.timer_start('main loop:readers')
        # Initialise ndarray to hold environment variables
        dtype = [(var, np.float32) for var in variables]
        env = np.ma.array(np.zeros(len(lon))*np.nan, dtype=dtype)

        if not hasattr(self, 'fallback_values'):
            self.set_fallback_values(refresh=False)

        # Discard any existing readers which are not relevant
        self.discard_irrelevant_readers()

        if 'drift:truncate_ocean_model_below_m' in self._config:
            truncate_depth = self.get_config('drift:truncate_ocean_model_below_m')
            if truncate_depth is not None:
                self.logger.debug('Truncating ocean models below %s m' % truncate_depth)
                z = z.copy()
                z[z<-truncate_depth] = -truncate_depth
                if self.required_profiles_z_range is not None:
                    self.required_profiles_z_range = np.array(
                        self.required_profiles_z_range)
                    self.required_profiles_z_range[self.required_profiles_z_range<-truncate_depth] = -truncate_depth

        # Initialise more lazy readers if necessary
        missing_variables = ['missingvar']
        while (len(missing_variables) > 0 and
               len(self._lazy_readers()) > 0):
            variable_groups, reader_groups, missing_variables = \
                self.get_reader_groups(variables)
            if hasattr(self, 'desired_variables'):
                missing_variables = list(set(missing_variables) -
                                         set(self.desired_variables))
            if len(missing_variables) > 0:
                self.logger.debug('Variables not covered by any reader: ' +
                              str(missing_variables))
                reader = 'NotNone'
                while reader is not None:
                    reader = self._initialise_next_lazy_reader()
                    if reader is not None:
                        if self.discard_reader_if_not_relevant(reader):
                            reader = None
                    if reader is not None:
                        if (reader.covers_time(self.time) and
                                len(reader.covers_positions(
                                lon, lat)[0]) > 0):
                            missing_variables = list(
                                set(missing_variables) -
                                set(reader.variables))
                            if len(missing_variables) == 0:
                                break  # We cover now all variables

        # For each variable/reader group:
        variable_groups, reader_groups, missing_variables = \
            self.get_reader_groups(variables)
        for variable in variables:  # Fill with fallback value if no reader
            co = self.get_config('environment:fallback:%s' % variable)
            if co is not None:
                env[variable] = np.ma.ones(env[variable].shape) * co

        for i, variable_group in enumerate(variable_groups):
            self.logger.debug('----------------------------------------')
            self.logger.debug('Variable group %s' % (str(variable_group)))
            self.logger.debug('----------------------------------------')
            reader_group = reader_groups[i]
            missing_indices = np.array(range(len(lon)))
            # For each reader:
            for reader_name in reader_group:
                self.logger.debug('Calling reader ' + reader_name)
                self.logger.debug('----------------------------------------')
                self.timer_start('main loop:readers:' +
                                 reader_name.replace(':', '<colon>'))
                reader = self.readers[reader_name]
                if reader.is_lazy:
                    self.logger.warning('Reader is lazy, should not happen')
                    import sys; sys.exit('Should not happen')
                if not reader.covers_time(time):
                    self.logger.debug('\tOutside time coverage of reader.')
                    if reader_name == reader_group[-1]:
                        if self._initialise_next_lazy_reader() is not None:
                            self.logger.debug('Missing variables: calling get_environment recursively')
                            return self.get_environment(variables,
                                        time, lon, lat, z, profiles)
                    continue
                # Fetch given variables at given positions from current reader
                try:
                    self.logger.debug('Data needed for %i elements' %
                                  len(missing_indices))
                    # Check if vertical profiles are requested from reader
                    if profiles is not None:
                        profiles_from_reader = list(
                            set(variable_group) & set(profiles))
                        if profiles_from_reader == []:
                            profiles_from_reader = None
                    else:
                        profiles_from_reader = None
                    env_tmp, env_profiles_tmp = \
                        reader.get_variables_interpolated(
                            variable_group, profiles_from_reader,
                            self.required_profiles_z_range, time,
                            lon[missing_indices], lat[missing_indices],
                            z[missing_indices], self.use_block, self.proj)

                except Exception as e:
                    self.logger.info('========================')
                    self.logger.info('Exception:')
                    self.logger.info(e)
                    self.logger.debug(traceback.format_exc())
                    self.logger.info('========================')
                    self.timer_end('main loop:readers:' +
                                   reader_name.replace(':', '<colon>'))
                    if reader_name == reader_group[-1]:
                        if self._initialise_next_lazy_reader() is not None:
                            self.logger.debug('Missing variables: calling get_environment recursively')
                            return self.get_environment(variables,
                                        time, lon, lat, z, profiles)
                    continue

                # Copy retrieved variables to env array, and mask nan-values
                for var in variable_group:
                    if var not in self.required_variables:
                        self.logger.debug('Not returning env-variable: ' + var)
                        continue
                    env[var][missing_indices] = np.ma.masked_invalid(
                        env_tmp[var][0:len(missing_indices)]).astype('float32')
                    if profiles_from_reader is not None and var in profiles_from_reader:
                        if 'env_profiles' not in locals():
                            env_profiles = env_profiles_tmp
                        # TODO: fix to be checked
                        if var in env_profiles and var in env_profiles_tmp:
                            # If one profile has fewer vertical layers than
                            # the other, we use only the overlapping part
                            if len(env_profiles['z']) != len(
                                env_profiles_tmp['z']):
                                self.logger.debug('Warning: different number of '
                                    ' vertical layers: %s and %s' % (
                                        len(env_profiles['z']),
                                        len( env_profiles_tmp['z'])))
                            z_ind = np.arange(np.minimum(
                                len(env_profiles['z'])-1,
                                len(env_profiles_tmp['z'])-1))
                            # len(missing_indices) since 2 points might have been added and not removed
                            env_profiles_tmp[var] = np.ma.atleast_2d(env_profiles_tmp[var])
                            env_profiles[var][np.ix_(z_ind, missing_indices)] = \
                                np.ma.masked_invalid(env_profiles_tmp[var][z_ind,0:len(missing_indices)]).astype('float32')
                            # For profiles with different numbers of layers, we extrapolate
                            if env_profiles[var].shape[0] > 1:
                                missingbottom = np.isnan(env_profiles[var][-1,:])
                                env_profiles[var][-1, missingbottom] = env_profiles[var][-2, missingbottom]

                # Detect elements with missing data, for present reader group
                if hasattr(env_tmp[variable_group[0]], 'mask'):
                    try:
                        del combined_mask
                    except:
                        pass
                    for var in variable_group:
                        tmp_var = np.ma.masked_invalid(env_tmp[var])
                        # Changed 13 Oct 2016, but uncertain of effect
                        # TODO: to be checked
                        #tmp_var = env_tmp[var]
                        if 'combined_mask' not in locals():
                            combined_mask = np.ma.getmask(tmp_var)
                        else:
                            combined_mask = \
                                np.ma.mask_or(combined_mask,
                                              np.ma.getmask(tmp_var),
                                              shrink=False)
                    try:
                        if len(missing_indices) != len(combined_mask):
                            # TODO: mask mismatch due to 2 added points
                            raise ValueError('Mismatch of masks')
                        missing_indices = missing_indices[combined_mask]
                    except:  # Not sure what is happening here
                        self.logger.info('Problems setting mask on missing_indices!')
                else:
                    missing_indices = []  # temporary workaround
                if (type(missing_indices) == np.int64) or (
                        type(missing_indices) == np.int32):
                    missing_indices = []
                self.timer_end('main loop:readers:' +
                               reader_name.replace(':', '<colon>'))
                if len(missing_indices) == 0:
                    self.logger.debug('Obtained data for all elements.')
                    break
                else:
                    self.logger.debug('Data missing for %i elements.' %
                                  (len(missing_indices)))
                    if len(self._lazy_readers()) > 0:
                        if self._initialise_next_lazy_reader() is not None:
                            self.logger.warning('Missing variables: calling get_environment recursively')
                            return self.get_environment(variables,
                                time, lon, lat, z, profiles)

        self.logger.debug('---------------------------------------')
        self.logger.debug('Finished processing all variable groups')

        self.timer_start('main loop:readers:postprocessing')
        for var in self.fallback_values:
            if (var not in variables) and (profiles is None or var not in profiles):
                continue
            mask = env[var].mask
            if any(mask==True):
                self.logger.debug('    Using fallback value %s for %s for %s elements' %
                              (self.fallback_values[var], var, np.sum(mask==True)))
                env[var][mask] = self.fallback_values[var]
            # Profiles
            if profiles is not None and var in profiles:
                if 'env_profiles' not in locals():
                    self.logger.debug('Creating empty dictionary for profiles not '
                                  'profided by any reader: ' + str(self.required_profiles))
                    env_profiles = {}
                    env_profiles['z'] = \
                        np.array(self.required_profiles_z_range)[::-1]
                if var not in env_profiles:
                    self.logger.debug('      Using fallback value %s for %s for all profiles' %
                                  (self.fallback_values[var], var))
                    env_profiles[var] = self.fallback_values[var]*\
                        np.ma.ones((len(env_profiles['z']), self.num_elements_active()))
                else:
                    mask = env_profiles[var].mask
                    num_masked_values_per_element = np.sum(mask==True)
                    num_missing_profiles = np.sum(num_masked_values_per_element == len(env_profiles['z']))
                    env_profiles[var][mask] = self.fallback_values[var]
                    self.logger.debug('      Using fallback value %s for %s for %s profiles' %
                                  (self.fallback_values[var], var, num_missing_profiles,))
                    num_missing_individual = np.sum(num_masked_values_per_element > 0) - num_missing_profiles
                    if num_missing_individual > 0:
                        self.logger.debug('        ...plus %s individual points in other profiles' %
                                      num_missing_individual)

        #######################################################
        # Some extra checks of units and realistic magnitude
        #######################################################
        if 'sea_water_temperature' in variables:
            t_kelvin = np.where(env['sea_water_temperature']>100)[0]
            if len(t_kelvin) > 0:
                self.logger.warning('Converting temperatures from Kelvin to Celcius')
                env['sea_water_temperature'][t_kelvin] = env['sea_water_temperature'][t_kelvin] - 273.15
                if 'env_profiles' in locals() and 'sea_water_temperature' in env_profiles.keys():
                  env_profiles['sea_water_temperature'][:,t_kelvin] = \
                    env_profiles['sea_water_temperature'][:,t_kelvin] - 273.15

        #######################################################
        # Parameterisation of unavailable variables
        #######################################################
        if 'drift:use_tabularised_stokes_drift' in self._config and self.get_config('drift:use_tabularised_stokes_drift') is True:
            if 'x_wind' not in variables:
                self.logger.debug('No wind available to calculate Stokes drift')
            else:
                if 'sea_surface_wave_stokes_drift_x_velocity' not in variables or (
                    env['sea_surface_wave_stokes_drift_x_velocity'].max() == 0 and
                    env['sea_surface_wave_stokes_drift_y_velocity'].max() == 0):
                        self.logger.debug('Calculating parameterised stokes drift')
                        env['sea_surface_wave_stokes_drift_x_velocity'], \
                        env['sea_surface_wave_stokes_drift_y_velocity'] = \
                            self.wave_stokes_drift_parameterised((env['x_wind'], env['y_wind']),
                                self.get_config('drift:tabularised_stokes_drift_fetch'))

                if (env['sea_surface_wave_significant_height'].max() == 0):
                        self.logger.debug('Calculating parameterised significant wave height')
                        env['sea_surface_wave_significant_height'] = \
                            self.wave_significant_height_parameterised((env['x_wind'], env['y_wind']),
                            self.get_config('drift:tabularised_stokes_drift_fetch'))

        #############################
        # Add uncertainty/diffusion
        #############################
        # Current
        if 'x_sea_water_velocity' in variables and \
                'y_sea_water_velocity' in variables:
            std = self.get_config('drift:current_uncertainty')
            if std > 0:
                self.logger.debug('Adding uncertainty for current: %s m/s' % std)
                env['x_sea_water_velocity'] += np.random.normal(
                    0, std, self.num_elements_active())
                env['y_sea_water_velocity'] += np.random.normal(
                    0, std, self.num_elements_active())
            std = self.get_config('drift:current_uncertainty_uniform')
            if std > 0:
                self.logger.debug('Adding uncertainty for current: %s m/s' % std)
                env['x_sea_water_velocity'] += np.random.uniform(
                    -std, std, self.num_elements_active())
                env['y_sea_water_velocity'] += np.random.uniform(
                    -std, std, self.num_elements_active())
        # Wind
        if 'x_wind' in variables and 'y_wind' in variables:
            std = self.get_config('drift:wind_uncertainty')
            if std > 0:
                self.logger.debug('Adding uncertainty for wind: %s m/s' % std)
                env['x_wind'] += np.random.normal(
                    0, std, self.num_elements_active())
                env['y_wind'] += np.random.normal(
                    0, std, self.num_elements_active())

        #####################
        # Diagnostic output
        #####################
        if len(env) > 0:
            self.logger.debug('------------ SUMMARY -------------')
            for var in variables:
                self.logger.debug('    %s: %g (min) %g (max)' %
                              (var, env[var].min(), env[var].max()))
            self.logger.debug('---------------------------------')
            self.logger.debug('\t\t%s active elements' % self.num_elements_active())
            if self.num_elements_active() > 0:
                lonmin = self.elements.lon.min()
                lonmax = self.elements.lon.max()
                latmin = self.elements.lat.min()
                latmax = self.elements.lat.max()
                zmin = self.elements.z.min()
                zmax = self.elements.z.max()
                if latmin == latmax:
                    self.logger.debug('\t\tlatitude =  %s' % (latmin))
                else:
                    self.logger.debug('\t\t%s <- latitude  -> %s' % (latmin, latmax))
                if lonmin == lonmax:
                    self.logger.debug('\t\tlongitude = %s' % (lonmin))
                else:
                    self.logger.debug('\t\t%s <- longitude -> %s' % (lonmin, lonmax))
                if zmin == zmax:
                    self.logger.debug('\t\tz = %s' % (zmin))
                else:
                    self.logger.debug('\t\t%s   <- z ->   %s' % (zmin, zmax))
                self.logger.debug('---------------------------------')

        # Prepare array indiciating which elements contain any invalid values
        missing = np.ma.masked_invalid(env[variables[0]]).mask
        for var in variables[1:]:
            missing = np.ma.mask_or(missing,
                                    np.ma.masked_invalid(env[var]).mask,
                                    shrink=False)

        # Convert dictionary to recarray and return
        if 'env_profiles' not in locals():
            env_profiles = None

        # Convert masked arrays to regular arrays for increased performance
        env = np.array(env)
        if env_profiles is not None:
            for var in env_profiles:
                env_profiles[var] = np.array(env_profiles[var])

        self.timer_end('main loop:readers:postprocessing')
        self.timer_end('main loop:readers')

        return env.view(np.recarray), env_profiles, missing

    def num_elements_active(self):
        """The number of active elements."""
        if hasattr(self, 'elements'):
            return len(self.elements)
        else:
            return 0

    def num_elements_deactivated(self):
        """The number of deactivated elements."""
        if hasattr(self, 'elements_deactivated'):
            return len(self.elements_deactivated)
        else:
            return 0

    def num_elements_scheduled(self):
        if hasattr(self, 'elements_scheduled'):
            return len(self.elements_scheduled)
        else:
            return 0

    def num_elements_total(self):
        """The total number of scheduled, active and deactivated elements."""
        return self.num_elements_activated() + self.num_elements_scheduled()

    def num_elements_activated(self):
        """The total number of active and deactivated elements."""
        return self.num_elements_active() + self.num_elements_deactivated()

    def schedule_elements(self, elements, time):
        """Schedule elements to be seeded during runtime.

        Also assigns a unique ID to each particle, monotonically increasing."""

        # prepare time
        if isinstance(time, np.ndarray):
            time = list(time)
        if not isinstance(time, list):
            time = [time]
        if len(time) == 1 and len(elements) > 1:
            time = time*len(elements)

        if not hasattr(self, 'elements_scheduled'):
            self.elements_scheduled = elements
            self.elements_scheduled_time = np.array(time)
            # We start simulation at time of release of first element:
            self.start_time = time[0]
            self.elements_scheduled.ID = np.arange(1, len(elements) + 1)
        else:
            elements.ID = np.arange(self.num_elements_scheduled() + 1,
                                    self.num_elements_scheduled() + 1 +
                                    len(elements))  # Increase ID successively
            self.elements_scheduled.extend(elements)
            self.elements_scheduled_time = np.append(
                self.elements_scheduled_time, np.array(time))

        min_time = np.min(time)
        if hasattr(self, 'start_time'):
            if min_time < self.start_time:
                self.start_time = min_time
                self.logger.debug('Setting simulation start time to %s' %
                              str(min_time))
        else:
            self.start_time = min_time
            self.logger.debug('Setting simulation start time to %s' %
                          str(min_time))

    def release_elements(self):
        """Activate elements which are scheduled within following timestep."""

        self.logger.debug('to be seeded: %s, already seeded %s' % (
            len(self.elements_scheduled), self.num_elements_activated()))
        if len(self.elements_scheduled) == 0:
            return
        if self.time_step.days >= 0:
            indices = (self.elements_scheduled_time >= self.time) & \
                      (self.elements_scheduled_time <
                       self.time + self.time_step)
        else:
            indices = (self.elements_scheduled_time <= self.time) & \
                      (self.elements_scheduled_time >
                       self.time + self.time_step)
        self.store_present_positions(
            self.elements_scheduled.ID[indices],
            self.elements_scheduled.lon[indices],
            self.elements_scheduled.lat[indices])
        self.elements_scheduled.move_elements(self.elements, indices)
        self.elements_scheduled_time = self.elements_scheduled_time[~indices]
        self.logger.debug('Released %i new elements.' % np.sum(indices))

    def closest_ocean_points(self, lon, lat):
        """Return the closest ocean points for given lon, lat"""

        deltalon = 0.01 # grid
        deltalat = 0.01
        numbuffer = 10
        lonmin = lon.min() - deltalon*numbuffer
        lonmax = lon.max() + deltalon*numbuffer
        latmin = lat.min() - deltalat*numbuffer
        latmax = lat.max() + deltalat*numbuffer
        if not 'land_binary_mask' in self.priority_list:
            self.logger.info('No land reader added, '
                         'making a temporary landmask reader')
            from opendrift.models.oceandrift import OceanDrift
            from opendrift.readers import reader_global_landmask
            reader_landmask = reader_global_landmask.Reader(
                    extent = [
                        np.maximum(-360, self.elements_scheduled.lon.min() - deltalon),
                        np.maximum(-89, self.elements_scheduled.lat.min() - deltalat),
                        np.minimum(720, self.elements_scheduled.lon.max() + deltalon),
                        np.minimum(89, self.elements_scheduled.lat.max() + deltalat)
                        ])
            reader_landmask.name = 'tempreader'
            o = OceanDrift(
                loglevel=self.logger.getEffectiveLevel())
            o.add_reader(reader_landmask)
            land_reader = reader_landmask

            tmp_reader = True
        else:
            self.logger.info('Using existing reader for land_binary_mask')
            land_reader_name = self.priority_list['land_binary_mask'][0]
            land_reader = self.readers[land_reader_name]
            o = self
            tmp_reader = False
        land = o.get_environment(['land_binary_mask'],
            lon=lon, lat=lat, z=0*lon, time=land_reader.start_time,
            profiles=None)[0]['land_binary_mask']
        if land.max() == 0:
            self.logger.info('All points are in ocean')
            return lon, lat
        self.logger.info('Moving %i out of %i points from land to water' %
                     (np.sum(land==1), len(lon)))
        landlons = lon[land==1]
        landlats = lat[land==1]
        longrid = np.arange(lonmin, lonmax, deltalon)
        latgrid = np.arange(latmin, latmax, deltalat)
        longrid, latgrid = np.meshgrid(longrid, latgrid)
        longrid = longrid.ravel()
        latgrid = latgrid.ravel()
        # Remove grid-points not covered by this reader
        latgrid_covered = land_reader.covers_positions(longrid, latgrid)[0]
        longrid = longrid[latgrid_covered]
        latgrid = latgrid[latgrid_covered]
        landgrid = o.get_environment(
            ['land_binary_mask'], lon=longrid, lat=latgrid,
            z=0*longrid, time=land_reader.start_time,
            profiles=None)[0]['land_binary_mask']
        if landgrid.min() == 1 or np.isnan(landgrid.min()):
            self.logger.warning('No ocean pixels nearby, cannot move elements.')
            return lon, lat

        oceangridlons = longrid[landgrid==0]
        oceangridlats = latgrid[landgrid==0]
        from scipy import spatial
        tree = scipy.spatial.cKDTree(
            np.dstack([oceangridlons, oceangridlats])[0])
        landpoints = np.dstack([landlons, landlats])
        dist, indices = tree.query(landpoints)
        indices = indices.ravel()
        lon[land==1] = oceangridlons[indices]
        lat[land==1] = oceangridlats[indices]

        return lon, lat

    def horizontal_diffusion(self):
        """Move elements with random walk according to given horizontal diffuivity."""
        D = self.get_config('drift:horizontal_diffusivity')
        if D == 0:
            self.logger.debug('Horizontal diffusivity is 0, no random walk.')
            return
        dt = self.time_step.total_seconds()
        x_vel = np.sqrt(2*D/dt)*np.random.normal(scale=1, size=self.num_elements_active())
        y_vel = np.sqrt(2*D/dt)*np.random.normal(scale=1, size=self.num_elements_active())
        speed = np.sqrt(x_vel*x_vel+y_vel*y_vel)
        self.logger.debug('Moving elements according to horizontal diffusivity of %s, with speeds between %s and %s m/s'
                          % (D, speed.min(), speed.max()))
        self.update_positions(x_vel, y_vel)

    def deactivate_elements(self, indices, reason='deactivated'):
        """Schedule deactivated particles for deletion (at end of step)"""
        if any(indices) is False:
            return
        if reason not in self.status_categories:
            self.status_categories.append(reason)
            self.logger.debug('Added status %s' % (reason))
        reason_number = self.status_categories.index(reason)
        #if not hasattr(self.elements.status, "__len__"):
        if len(np.atleast_1d(self.elements.status)) == 1:
            status = self.elements.status.item()
            self.elements.status = np.zeros(self.num_elements_active())
            self.elements.status.fill(status)
        # Deactivate elements, if they have not already been deactivated
        self.elements.status[indices & (self.elements.status ==0)] = \
            reason_number
        self.logger.debug('%s elements scheduled for deactivation (%s)' %
                      (np.sum(indices), reason))
        self.logger.debug('\t(z: %f to %f)' %
            (self.elements.z[indices].min(),
             self.elements.z[indices].max()))

    def remove_deactivated_elements(self):
        """Moving deactivated elements from self.elements
        to self.elements_deactivated."""

        # All particles scheduled for deletion
        indices = (self.elements.status != 0)
        #try:
        #    len(indices)
        #except:
        if len(indices) == 0 or np.sum(indices) == 0:
            self.logger.debug('No elements to deactivate')
            return  # No elements scheduled for deactivation
        # Basic, but some more housekeeping will be required later
        self.elements.move_elements(self.elements_deactivated, indices)
        self.logger.debug('Removed %i elements.' % (np.sum(indices)))
        if hasattr(self, 'environment'):
            self.environment = self.environment[~indices]
            self.logger.debug('Removed %i values from environment.' %
                          (np.sum(indices)))
        if hasattr(self, 'environment_profiles') and \
                self.environment_profiles is not None:
            for varname, profiles in self.environment_profiles.items():
                self.logger.debug('remove items from profile for '+varname)
                if varname != 'z':
                    self.environment_profiles[varname] = \
                        profiles[:, ~indices]
            self.logger.debug('Removed %i values from environment_profiles.' %
                          (np.sum(indices)))
            #if self.num_elements_active() == 0:
            #    raise ValueError('No more active elements.')  # End simulation

    def set_fallback_values(self, refresh=False):
        if hasattr(self, 'fallback_values') and refresh is False:
            raise ValueError('Manually editing fallback_values dict is deprecated, please use set_config()')
        else:
            c = self.get_configspec('environment:fallback:')
            self.fallback_values = {}
            for var in list(c):
                if c[var]['value'] is not None:
                    self.fallback_values[var.split(':')[-1]] = c[var]['value']

    def run(self, time_step=None, steps=None, time_step_output=None,
            duration=None, end_time=None, outfile=None, export_variables=None,
            export_buffer_length=100, stop_on_error=False):
        """Start a trajectory simulation, after initial configuration.

        Performs the main loop:
            - Obtain environment data for positions of all particles.
            - Call method 'update' to update (incl advect) particle properties.
        until one of the following conditions are met:
            - Maximum number of steps are reached
            - A needed variable can not be obtained by any reader
                (outside spatial/temporal domain) and has no fallback
                (default) value.
            - All particles have been deactivated (e.g. by stranding)
            - Occurance of any error, whose trace will be output to terminal.

        Before starting a model run, readers must be added for all
        required variables, unless fallback values have been specified.
        Some particles/elements must have been scheduled for seeding, and the
        run will start at the time when the first element has been scheduled..

        Arguments:
            time_step: interval between particles updates, in seconds or as
                timedelta. Default: 3600 seconds (1 hour)
            time_step_output: Time step at which element properties are stored
                and eventually written to file.
                Timedelta object or seconds.
                Default: same as time_step, meaning that all steps are stored
            The length of the simulation is specified by defining one
                (and only one) of the following parameters:
                - steps: integer, maximum number of steps. End of simulation
                will be self.start_time + steps*self.time_step
                - duration: timedelta defining the length of the simulation
                - end_time: datetime object defining the end of the simulation
            export_variables: list of variables and parameter names to be
                saved to file. Default is None (all variables are saved)
        """

        # Exporting software and hardware specification, for possible debugging
        self.logger.debug(opendrift.versions())

        self.timer_end('configuration')
        self.timer_start('preparing main loop')

        if self.num_elements_scheduled() == 0:
            raise ValueError('Please seed elements before starting a run.')
        self.elements = self.ElementType()

        # Collect fallback values from config into dict
        self.set_fallback_values(refresh=True)

        if outfile is None and export_buffer_length is not None:
            self.logger.debug('No output file is specified, '
                          'neglecting export_buffer_length')
            export_buffer_length = None

        # Set projection to latlong if not taken from any of the readers
        if self.proj is not None and not (self.proj.crs.is_geographic or
            'proj=merc' in self.proj.srs):
            for vector_component in vector_pairs_xy:
                for component in vector_component:
                    if component in self.fallback_values and \
                            self.fallback_values[component] != 0:
                        self.logger.info('Setting SRS to latlong, since non-zero '
                                     'value used for fallback vectors (%s)' %
                                     component)
                        self.set_projection('+proj=latlong')
        if self.proj is None:
            self.logger.info('Setting SRS to latlong, since not defined before.')
            self.set_projection('+proj=latlong')

        # Check if any readers have same SRS as simulation
        for reader in self.readers.values():
            if reader.is_lazy:
                continue
            readerSRS = reader.proj.srs.replace(' +ellps=WGS84', '').strip()
            simulationSRS = self.proj.srs.replace(' +ellps=WGS84', '').strip()
            if readerSRS == simulationSRS:
                reader.simulation_SRS = True
            else:
                reader.simulation_SRS = False

        # Make constant readers if config environment:constant:<var> is
        c = self.get_configspec('environment:constant:')
        mr = {}
        for var in list(c):
            if c[var]['value'] is not None:
                mr[var.split(':')[-1]] = c[var]['value']
        if len(mr) > 0:
            from opendrift.readers import reader_constant
            rc = reader_constant.Reader(mr)
            self.add_reader(rc, first=True)

        missing_variables = self.missing_variables()
        missing_variables = [m for m in missing_variables if
                             m != 'land_binary_mask']
        if len(missing_variables) > 0:
            has_fallback = [var for var in missing_variables
                            if var in self.fallback_values]
            has_no_fallback = [var for var in missing_variables
                               if var not in self.fallback_values]
            #if has_fallback == missing_variables:
            if len(has_fallback) > 0:# == missing_variables:
                self.logger.info('Fallback values will be used for the following '
                             'variables which have no readers: ')
                for var in has_fallback:
                    self.logger.info('\t%s: %f' % (var, self.fallback_values[var]))
            #else:
            if len(has_no_fallback) > 0 and len(self._lazy_readers()) == 0:# == missing_variables:
                self.logger.warning('No readers added for the following variables: '
                                + str(has_no_fallback))
                raise ValueError('Readers must be added for the '
                                 'following required variables: ' +
                                 str(has_no_fallback))

        # Some cleanup needed if starting from imported state
        if self.steps_calculation >= 1:
            self.steps_calculation = 0
        if hasattr(self, 'history'):
            # Delete history matrix before new run
            delattr(self, 'history')
            # Renumbering elements from 0 to num_elements, necessary fix when
            # importing from file, where elements may have been deactivated
            # TODO: should start from 1?
            self.elements.ID = np.arange(0, self.num_elements_active())

        ########################
        # Simulation time step
        ########################
        if time_step is None:
            time_step = timedelta(minutes=self.get_config('general:time_step_minutes'))
        if type(time_step) is not timedelta:
            # Time step may be given in seconds, as alternative to timedelta
            time_step = timedelta(seconds=time_step)
        self.time_step = time_step
        if time_step_output is None:
            time_step_output = self.get_config('general:time_step_output_minutes')
            if time_step_output is None:
                self.time_step_output = self.time_step
            else:
                self.time_step_output = timedelta(minutes=time_step_output)
        else:
            if type(time_step_output) is timedelta:
                self.time_step_output = time_step_output
            else:
                self.time_step_output = timedelta(seconds=time_step_output)
            if self.time_step_output.days >= 0 and self.time_step.days < 0:
                self.time_step_output = -self.time_step_output

        time_step_ratio = self.time_step_output.total_seconds() / \
            self.time_step.total_seconds()
        if time_step_ratio < 1:
            raise ValueError('Output time step must be equal or larger '
                             'than calculation time step.')
        if not time_step_ratio.is_integer():
            raise ValueError('Ratio of calculation and output time steps '
                             'must be an integer - given ratio is %s' %
                             time_step_ratio)
        ########################
        # Simulation duration
        ########################
        if time_step.days < 0:
            self.logger.info('Backwards simulation, starting from last seeded element')
            self.start_time = self.elements_scheduled_time.max()
        if (duration is not None and end_time is not None) or \
            (duration is not None and steps is not None) or \
                (steps is not None and end_time is not None):
            raise ValueError('Only one of "steps", "duration" and "end_time" '
                             'may be provided simultaneously')
        if duration is None and end_time is None:
            if steps is not None:
                duration = steps*self.time_step
            else:
                for reader in self.readers.values():
                    if reader.end_time is not None:
                        if end_time is None:
                            end_time = reader.end_time
                        else:
                            end_time = min(end_time, reader.end_time)
                    self.logger.info('Duration, steps or end time not specified, '
                                 'running until end of first reader: %s' %
                                 (end_time))
        if duration is None:
            duration = end_time - self.start_time

        if time_step.days < 0 and duration.days >= 0:
            # Duration shall also be negative for backwards run
            duration = -duration

        if np.sign(duration.total_seconds()) * np.sign(time_step.total_seconds()) < 0:
            raise ValueError("Time step must be negative if duration is negative.")

        self.expected_steps_output = duration.total_seconds() / \
            self.time_step_output.total_seconds() + 1  # Includes start and end
        self.expected_steps_calculation = duration.total_seconds() / \
            self.time_step.total_seconds()
        self.expected_steps_output = int(self.expected_steps_output)
        self.expected_steps_calculation = int(self.expected_steps_calculation)
        self.expected_end_time = self.start_time + self.expected_steps_calculation*self.time_step

        ##############################################################
        # Prepare readers for the requested simulation domain/time
        ##############################################################
        max_distance = \
            self.max_speed*self.expected_steps_calculation * \
            np.abs(self.time_step.total_seconds())
        deltalat = max_distance/111000.
        deltalon = deltalat/np.cos(
            np.radians(np.mean(self.elements_scheduled.lat)))
        # TODO: extent should ideally be a general polygon, not only lon/lat-min/max
        # TODO: Should also take into account eventual lifetime of elements
        simulation_extent = [
            np.maximum(-360, self.elements_scheduled.lon.min() - deltalon),
            np.maximum(-89, self.elements_scheduled.lat.min() - deltalat),
            np.minimum(720, self.elements_scheduled.lon.max() + deltalon),
            np.minimum(89, self.elements_scheduled.lat.max() + deltalat)]
        self.logger.debug('Preparing readers for simulation coverage (%s) and time (%s to %s)'
                % (simulation_extent, self.start_time, self.expected_end_time))
        for reader in self.readers.values():
            self.logger.debug('\tPreparing %s' % reader.name)
            reader.prepare(
                extent=simulation_extent,
                start_time=self.start_time, end_time = self.expected_end_time)

        ##############################################################
        # If no landmask has been added, we determine it dynamically
        ##############################################################
        # TODO: some more error checking here
        # If landmask is requested, it shall not be obtained from other readers
        if self.get_config('general:use_auto_landmask') is True:
            if 'land_binary_mask' in self.priority_list:
                if 'basemap_landmask' in self.priority_list['land_binary_mask']:
                    self.priority_list['land_binary_mask'] = ['basemap_landmask']
                elif 'global_landmask' in self.priority_list['land_binary_mask']:
                    self.priority_list['land_binary_mask'] = ['global_landmask']
                else:
                    del self.priority_list['land_binary_mask']
        if self.get_config('general:use_auto_landmask') is True and \
                ('land_binary_mask' in self.required_variables and \
                'land_binary_mask' not in self.priority_list \
                and 'land_binary_mask' not in self.fallback_values):
            self.logger.info(
                'Adding a dynamical landmask with max. priority based on '
                'assumed maximum speed of %s m/s. '
                'Adding a customised landmask may be faster...' % self.max_speed)
            self.timer_start('preparing main loop:making dynamical landmask')
            from opendrift.readers import reader_global_landmask
            reader_landmask = reader_global_landmask.Reader(extent = simulation_extent)
            self.add_reader(reader_landmask)

            self.timer_end('preparing main loop:making dynamical landmask')

        # Move point seed on land to ocean
        if self.get_config('seed:ocean_only') is True and \
            ('land_binary_mask' not in self.fallback_values) and \
            ('land_binary_mask' in self.required_variables):
            self.timer_start('preparing main loop:moving elements to ocean')
            self.elements_scheduled.lon, self.elements_scheduled.lat = \
                self.closest_ocean_points(self.elements_scheduled.lon,
                                          self.elements_scheduled.lat)
            self.timer_end('preparing main loop:moving elements to ocean')

        ####################################################################
        # Preparing history array for storage in memory and eventually file
        ####################################################################
        if export_buffer_length is None:
            self.export_buffer_length = self.expected_steps_output
        else:
            self.export_buffer_length = export_buffer_length

        if self.time_step.days < 0:
            # For backwards simulation, we start at last seeded element
            self.logger.info('Backwards simulation, starting at '
                         'time of last seeded element')
            self.time = self.elements_scheduled_time.max()
            # Flipping ID array, so that lowest IDs are released first
            self.elements_scheduled.ID = \
                np.flipud(self.elements_scheduled.ID)
        else:
            # Forward simulation, start time has been set when seeding
            self.time = self.start_time

        # Add the output variables which are always required
        if export_variables is not None:
            export_variables = list(set(export_variables +
                                        ['lon', 'lat', 'ID', 'status']))
        self.export_variables = export_variables
        # Initialise array to hold history (element properties and environment)
        # for export to file.
        history_dtype_fields = [(name,
                                 self.ElementType.variables[name]['dtype'])
                                for name in self.ElementType.variables]
        # Add environment variables
        self.history_metadata = self.ElementType.variables.copy()
        for env_var in self.required_variables:
            history_dtype_fields.append((env_var, np.dtype('float32')))
            self.history_metadata[env_var] = {}

        # Remove variables from output array, if only subset is requested
        if self.export_variables is not None:
            history_dtype_fields = [f for f in history_dtype_fields
                                    if f[0] in self.export_variables]
            for m in list(self.history_metadata):
                if m not in self.export_variables:
                    del self.history_metadata[m]

        history_dtype = np.dtype(history_dtype_fields)
        self.history = np.ma.array(np.zeros((len(self.elements_scheduled),
                                             self.export_buffer_length)),
                                   dtype=history_dtype)
        self.history.mask = True
        self.steps_exported = 0

        if outfile is not None:
            self.io_init(outfile)
        else:
            self.outfile = None

        #############################
        # Check validity domain
        #############################
        validity_domain = [
            self.get_config('drift:deactivate_west_of'),
            self.get_config('drift:deactivate_east_of'),
            self.get_config('drift:deactivate_south_of'),
            self.get_config('drift:deactivate_north_of')]
        if validity_domain == [None, None, None, None]:
            self.validity_domain = None
        else:
            self.validity_domain = validity_domain

        #############################
        # Model specific preparation
        #############################
        self.prepare_run()

        ##########################
        # Main loop
        ##########################
        self.add_metadata('simulation_time', datetime.now())
        self.timer_end('preparing main loop')
        self.timer_start('main loop')
        for i in range(self.expected_steps_calculation):
            try:
                # Release elements
                self.release_elements()

                if self.num_elements_active() == 0 and self.num_elements_scheduled() > 0:
                    self.steps_calculation += 1
                    self.logger.info('No active but %s scheduled elements, skipping timestep %s (%s)'
                                     % (self.num_elements_scheduled(), self.steps_calculation, self.time))
                    self.state_to_buffer()  # Append status to history array
                    if self.time is not None:
                        self.time = self.time + self.time_step
                    continue

                self.increase_age_and_retire()

                self.lift_elements_to_seafloor()  # If seafloor is penetrated

                if self.show_continuous_performance is True:
                    self.logger.info(self.performance())
                # Display time to terminal
                self.logger.debug('==================================='*2)
                self.logger.info('%s - step %i of %i - %i active elements '
                             '(%i deactivated)' %
                             (self.time, self.steps_calculation + 1,
                              self.expected_steps_calculation,
                              self.num_elements_active(),
                              self.num_elements_deactivated()))
                self.logger.debug('%s elements scheduled.' %
                              self.num_elements_scheduled())
                self.logger.debug('==================================='*2)

                self.environment, self.environment_profiles, missing = \
                    self.get_environment(list(self.required_variables),
                                         self.time,
                                         self.elements.lon,
                                         self.elements.lat,
                                         self.elements.z,
                                         self.required_profiles)

                self.store_previous_variables()

                self.calculate_missing_environment_variables()

                if any(missing):
                    self.report_missing_variables()

                self.interact_with_coastline()

                self.lift_elements_to_seafloor()  # If seafloor is penetrated

                self.deactivate_elements(missing, reason='missing_data')

                self.state_to_buffer()  # Append status to history array

                self.remove_deactivated_elements()

                # Propagate one timestep forwards
                self.steps_calculation += 1

                if self.num_elements_active() == 0 and self.num_elements_scheduled() == 0:
                    raise ValueError('No more active or scheduled elements, quitting.')

                # Store location, in case elements shall be moved back
                self.store_present_positions()

                #####################################################
                if self.num_elements_active() > 0:
                    self.logger.debug('Calling %s.update()' %
                                  type(self).__name__)
                    self.timer_start('main loop:updating elements')
                    self.update()
                    self.timer_end('main loop:updating elements')
                else:
                    self.logger.info('No active elements, skipping update() method')
                #####################################################

                self.horizontal_diffusion()

                if self.num_elements_active() == 0 and self.num_elements_scheduled() == 0:
                    raise ValueError('No active or scheduled elements, quitting simulation')

                self.logger.debug('%s active elements (%s deactivated)' %
                              (self.num_elements_active(),
                               self.num_elements_deactivated()))
                # Updating time
                if self.time is not None:
                    self.time = self.time + self.time_step

            except Exception as e:
                message = ('The simulation stopped before requested '
                           'end time was reached.')
                self.logger.warning(message)
                self.store_message(message)
                self.logger.info('========================')
                self.logger.info('End of simulation:')
                self.logger.info(e)
                self.logger.info(traceback.format_exc())
                self.logger.info(self.get_messages())
                if not hasattr(self, 'environment'):
                    sys.exit('Simulation aborted. ' +
                             self.get_messages())
                self.logger.info('========================')
                if stop_on_error is True:
                    sys.exit('Stopping on error. ' +
                             self.get_messages())
                if self.steps_calculation <= 1:
                    raise ValueError('Simulation stopped within '
                        'first timestep. ' + self.get_messages())
                break

        self.timer_end('main loop')
        self.timer_start('cleaning up')
        self.logger.debug('Cleaning up')

        self.interact_with_coastline(final=True)
        self.state_to_buffer()  # Append final status to buffer

        #############################
        # Add some metadata
        #############################
        for var in self.required_variables:
            keyword = 'reader_' + var
            if var not in self.priority_list:
                self.add_metadata(keyword, self.fallback_values[var])
            else:
                readers = self.priority_list[var]
                if readers[0].startswith('constant_reader'):
                    self.add_metadata(keyword, self.readers[readers[
                                0]]._parameter_value_map[var][0])
                else:
                    self.add_metadata(keyword,
                                      self.priority_list[var])

        if outfile is not None:
            self.logger.debug('Writing and closing output file: %s' % outfile)
            # Write buffer to outfile, and close
            if self.steps_output >= self.steps_exported:
                # Write last lines, if needed
                self.io_write_buffer()
            self.io_close()

        # Remove any elements scheduled for deactivation during last step
        #self.remove_deactivated_elements()

        if export_buffer_length is None:
            # Remove columns for unseeded elements in history array
            if self.num_elements_scheduled() > 0:
                self.logger.info('Removing %i unseeded elements from history array' %
                               self.num_elements_scheduled())
                mask = np.ones(self.history.shape[0], dtype=bool)
                mask[self.elements_scheduled.ID-1] = False
                self.history = self.history[mask, :]

            # Remove rows for unreached timsteps in history array
            self.history = self.history[:, range(self.steps_output)]
        else:  # If output has been flushed to file during run, we
               # need to reimport from file to get all data in memory
            del self.environment
            if hasattr(self, 'environment_profiles'):
                del self.environment_profiles
            self.io_import_file(outfile)

        self.timer_end('cleaning up')
        self.timer_end('total time')

    def increase_age_and_retire(self):
        """Increase age of elements, and retire if older than config setting."""
        # Increase age of elements
        self.elements.age_seconds += self.time_step.total_seconds()

        # Deactivate elements that exceed a certain age
        if self.get_config('drift:max_age_seconds') is not None:
            self.deactivate_elements(self.elements.age_seconds >=
                                     self.get_config('drift:max_age_seconds'),
                                     reason='retired')

        # Deacticate any elements outside validity domain set by user
        if self.validity_domain is not None:
            W, E, S, N = self.validity_domain
            if W is not None:
                self.deactivate_elements(self.elements.lon < W, reason='outside')
            if E is not None:
                self.deactivate_elements(self.elements.lon > E, reason='outside')
            if S is not None:
                self.deactivate_elements(self.elements.lat < S, reason='outside')
            if N is not None:
                self.deactivate_elements(self.elements.lat > N, reason='outside')

    def state_to_buffer(self):
        """Append present state (elements and environment) to recarray."""

        steps_calculation_float = \
            (self.steps_calculation * self.time_step.total_seconds() /
             self.time_step_output.total_seconds()) + 1
        if self.time_step <= timedelta(seconds=1):
            self.steps_output = int(np.round(steps_calculation_float))
        else:
            self.steps_output = int(np.floor(steps_calculation_float))

        ID_ind = self.elements.ID - 1
        time_ind = self.steps_output - 1 - self.steps_exported

        if steps_calculation_float.is_integer() or self.time_step < timedelta(seconds=1):
            element_ind = range(len(ID_ind))  # We write all elements
        else:
            deactivated = np.where(self.elements.status != 0)[0]
            if len(deactivated) == 0:
                    return  # No deactivated elements this sub-timestep
            # We write history for deactivated elements only:
            self.logger.debug('Writing history for %s deactivated elements' %
                          len(deactivated))
            ID_ind = ID_ind[deactivated]
            element_ind = deactivated
            time_ind = np.minimum(time_ind + 1, self.history.shape[1] - 1)

        # Store present state in history recarray
        for i, var in enumerate(self.elements.variables):
            if self.export_variables is not None and \
                    var not in self.export_variables:
                continue
            # Temporarily assuming elements numbered
            # from 0 to num_elements_active()
            # Does not hold when importing ID from a saved file, where
            # some elements have been deactivated
            self.history[var][ID_ind, time_ind] = \
                getattr(self.elements, var)[element_ind]
        # Copy environment data to history array
        for i, var in enumerate(self.environment.dtype.names):
            if self.export_variables is not None and \
                    var not in self.export_variables:
                continue
            self.history[var][ID_ind, time_ind] = \
                getattr(self.environment, var)[element_ind]

        # Call writer if buffer is full
        if (self.outfile is not None) and \
                ((self.steps_output - self.steps_exported) ==
                    self.export_buffer_length):
            self.io_write_buffer()

    def report_missing_variables(self):
        """Issue warning if some environment variables missing."""

        missing_variables = []
        for var in self.required_variables:
            if np.isnan(getattr(self.environment, var).min()):
                missing_variables.append(var)

        if len(missing_variables) > 0:
            self.logger.warning('Missing variables: ' +
                            str(missing_variables))
            self.store_message('Missing variables: ' +
                               str(missing_variables))

    def index_of_activation_and_deactivation(self):
        """Return the indices when elements were seeded and deactivated."""

        firstlast = np.ma.notmasked_edges(self.history['lon'], axis=1)
        index_of_activation = firstlast[0][1]
        index_of_deactivation = firstlast[1][1]
        if len(index_of_deactivation) < self.history['lon'].shape[0]:
            missingind = np.setdiff1d(
                np.arange(0, self.history['lon'].shape[0]),
                firstlast[0][0])
            self.logger.warning('%s elements were never seeded, removing from history array (this is probably caused by importing an old file)' % len(missingind))
            self.history = self.history[firstlast[0][0], :]

        return index_of_activation, index_of_deactivation

    def get_lonlats(self):
        if hasattr(self, 'history'):
            lons = self.history['lon']
            lats = self.history['lat']
        else:
            if self.steps_output > 0:
                lons = np.ma.array(np.reshape(self.elements.lon, (1, -1))).T
                lats = np.ma.array(np.reshape(self.elements.lat, (1, -1))).T
            else:
                lons = np.ma.array(
                    np.reshape(self.elements_scheduled.lon, (1, -1))).T
                lats = np.ma.array(
                    np.reshape(self.elements_scheduled.lat, (1, -1))).T
        return lons, lats

    def get_density_xarray(self, pixelsize_m, weight=None, per_origin_marker=False):
        from xhistogram.xarray import histogram

        if os.path.exists(self.analysis_file) and 'lon_bin' in xr.open_dataset(self.analysis_file):  # Import from file
            self.logger.info('Importing from saved analysis file')
            self.af = xr.open_dataset(self.analysis_file)
            histograms = []
            for var in self.af.variables:
                if var[0:9] == 'histogram':
                    histograms.append(self.af[var])
            lonbin = self.af.lon_bin.values
            latbin = self.af.lat_bin.values
            deltalon = (lonbin[1]-lonbin[0])
            lonbin = np.append(lonbin, lonbin[-1] + deltalon) - deltalon/2
            deltalat = (latbin[1]-latbin[0])
            latbin = np.append(latbin, latbin[-1] + deltalat) - deltalat/2
            self.af.close()

        else:  # Calculate
            start = datetime.now()
            self.logger.info('Calculating density array, this may take some time...')
            deltalat = pixelsize_m/111000.0  # m to degrees
            if hasattr(self, 'lonmin'):
                lonmin = self.lonmin
                lonmax = self.lonmax
                latmin = self.latmin
                latmax = self.latmax
            else:
                self.logger.debug('Finding min and max of lon and lat...')
                lonmin = np.nanmin(self.ds.lon)
                lonmax = np.nanmax(self.ds.lon)
                latmin = np.nanmin(self.ds.lat)
                latmax = np.nanmax(self.ds.lat)
                if os.path.exists(self.analysis_file):
                    self.af = Dataset(self.analysis_file, 'a')
                else:
                    self.af = Dataset(self.analysis_file, 'w')
                self.af.lonmin = lonmin
                self.af.lonmax = lonmax
                self.af.latmin = latmin
                self.af.latmax = latmax
                self.af.close()

            deltalon = deltalat/np.cos(np.radians((latmin+latmax)/2))
            latbin = np.arange(latmin-deltalat, latmax+deltalat, deltalat)
            lonbin = np.arange(lonmin-deltalon, lonmax+deltalon, deltalon)

            histograms = []
            if per_origin_marker is True:
                max_om = self.ds.origin_marker.max().compute().values
                for om in range(max_om+1):
                    self.logger.info('\tcalculating for origin_marker %s...' % om)
                    h = histogram(self.ds.lon.where(self.ds.origin_marker==om),
                                  self.ds.lat.where(self.ds.origin_marker==om), bins=[lonbin, latbin],
                                  dim=['trajectory'])
                    h.name = 'histogram_lon_lat_' + str(om)
                    histograms.append(h)
            else:
                h = histogram(self.ds.lon, self.ds.lat, bins=[lonbin, latbin],
                              dim=['trajectory'])
                histograms = [h]

            if self.analysis_file is not None:
                self.logger.info('Writing density array to analysis file: %s'
                                 % self.analysis_file)
                for om, h in enumerate(histograms):
                    if os.path.exists(self.analysis_file):
                        mode = 'a'
                    else:
                        mode = 'w'
                    self.logger.info('\twriting for origin_marker %s...' % om)
                    h.to_netcdf(self.analysis_file, mode)
            self.logger.info('Time to calculate density array: %s' %
                             (datetime.now() - start))

        return histograms, lonbin, latbin

    def get_density_array(self, pixelsize_m, weight=None):
        lon = self.get_property('lon')[0]
        lat = self.get_property('lat')[0]
        times = self.get_time_array()[0]
        deltalat = pixelsize_m/111000.0  # m to degrees
        deltalon = deltalat/np.cos(np.radians((np.nanmin(lat) +
                                               np.nanmax(lat))/2))
        lat_array = np.arange(np.nanmin(lat)-deltalat,
                              np.nanmax(lat)+deltalat, deltalat)
        lon_array = np.arange(np.nanmin(lon)-deltalat,
                              np.nanmax(lon)+deltalon, deltalon)
        bins=(lon_array, lat_array)
        z = self.get_property('z')[0]
        if weight is not None:
            weight_array = self.get_property(weight)[0]

        status = self.get_property('status')[0]
        lon_submerged = lon.copy()
        lat_submerged = lat.copy()
        lon_stranded = lon.copy()
        lat_stranded = lat.copy()
        lon_submerged[z>=0] = 1000
        lat_submerged[z>=0] = 1000
        lon[z<0] = 1000
        lat[z<0] = 1000
        H = np.zeros((len(times), len(lon_array) - 1,
                      len(lat_array) - 1))#.astype(int)
        H_submerged = H.copy()
        H_stranded = H.copy()
        try:
            strandnum = self.status_categories.index('stranded')
            lon_stranded[status!=strandnum] = 1000
            lat_stranded[status!=strandnum] = 1000
            contains_stranded = True
        except ValueError:
            contains_stranded = False

        for i in range(len(times)):
            if weight is not None:
                weights = weight_array[i,:]
            else:
                weights = None
            H[i,:,:], dummy, dummy = \
                np.histogram2d(lon[i,:], lat[i,:],
                               weights=weights, bins=bins)
            H_submerged[i,:,:], dummy, dummy = \
                np.histogram2d(lon_submerged[i,:], lat_submerged[i,:],
                               weights=weights, bins=bins)
            if contains_stranded is True:
                H_stranded[i,:,:], dummy, dummy = \
                np.histogram2d(lon_stranded[i,:], lat_stranded[i,:],
                               weights=weights, bins=bins)

        return H, H_submerged, H_stranded, lon_array, lat_array

    def get_residence_time(self, pixelsize_m):
        H,H_sub, H_str,lon_array,lat_array = \
            self.get_density_array(pixelsize_m)
        residence = np.sum(H, axis=0)
        return residence, lon_array, lat_array

    def write_netcdf_density_map(self, filename, pixelsize_m='auto'):
        '''Write netCDF file with map of particles densities'''

        if pixelsize_m == 'auto':
            lon, lat = self.get_lonlats()
            latspan = lat.max()-lat.min()
            pixelsize_m=30
            if latspan > .05:
                pixelsize_m = 50
            if latspan > .1:
                pixelsize_m = 300
            if latspan > .3:
                pixelsize_m = 500
            if latspan > .7:
                pixelsize_m = 1000
            if latspan > 2:
                pixelsize_m = 2000
            if latspan > 5:
                pixelsize_m = 4000

        H, H_submerged, H_stranded, lon_array, lat_array = \
            self.get_density_array(pixelsize_m)
        lon_array = (lon_array[0:-1] + lon_array[1::])/2
        lat_array = (lat_array[0:-1] + lat_array[1::])/2

        nc = Dataset(filename, 'w')
        nc.createDimension('lon', len(lon_array))
        nc.createDimension('lat', len(lat_array))
        nc.createDimension('time', H.shape[0])
        times = self.get_time_array()[0]
        timestr = 'seconds since 1970-01-01 00:00:00'
        nc.createVariable('time', 'f8', ('time',))
        nc.variables['time'][:] = date2num(times, timestr)
        nc.variables['time'].units = timestr
        nc.variables['time'].standard_name = 'time'
        # Projection
        nc.createVariable('projection_lonlat', 'i8')
        nc.variables['projection_lonlat'].grid_mapping_name = \
            'latitude_longitude'
        nc.variables['projection_lonlat'].earth_radius = 6371229.
        nc.variables['projection_lonlat'].proj4 = \
            '+proj=longlat +a=6371229 +no_defs'
        # Coordinates
        nc.createVariable('lon', 'f8', ('lon',))
        nc.createVariable('lat', 'f8', ('lat',))
        nc.variables['lon'][:] = lon_array
        nc.variables['lon'].long_name = 'longitude'
        nc.variables['lon'].short_name = 'longitude'
        nc.variables['lon'].units = 'degrees_east'
        nc.variables['lat'][:] = lat_array
        nc.variables['lat'].long_name = 'latitude'
        nc.variables['lat'].short_name = 'latitude'
        nc.variables['lat'].units = 'degrees_north'
        # Density
        nc.createVariable('density_surface', 'u1',
                          ('time','lat', 'lon'))
        H = np.swapaxes(H, 1, 2).astype('uint8')
        H = np.ma.masked_where(H==0, H)
        nc.variables['density_surface'][:] = H
        nc.variables['density_surface'].long_name = 'Detection probability'
        nc.variables['density_surface'].grid_mapping = 'projection_lonlat'
        nc.variables['density_surface'].units = '1'
        # Density submerged
        nc.createVariable('density_submerged', 'u1',
                          ('time','lat', 'lon'))
        H_sub = np.swapaxes(H_submerged, 1, 2).astype('uint8')
        H_sub = np.ma.masked_where(H_sub==0, H_sub)
        nc.variables['density_submerged'][:] = H_sub
        nc.variables['density_submerged'].long_name = 'Detection probability submerged'
        nc.variables['density_submerged'].grid_mapping = 'projection_lonlat'
        nc.variables['density_submerged'].units = '1'
        # Density stranded
        nc.createVariable('density_stranded', 'u1',
                          ('time','lat', 'lon'))
        H_stranded = np.swapaxes(H_stranded, 1, 2).astype('uint8')
        H_stranded = np.ma.masked_where(H_stranded==0, H_stranded)
        nc.variables['density_stranded'][:] = H_stranded
        nc.variables['density_stranded'].long_name = 'Detection probability stranded'
        nc.variables['density_stranded'].grid_mapping = 'projection_lonlat'
        nc.variables['density_stranded'].units = '1'

        nc.close()

    def get_time_array(self):
        """Return a list of output times of last run."""

        # Making sure start_time is datetime, and not cftime object
        self.start_time = datetime(self.start_time.year, self.start_time.month,
                                   self.start_time.day, self.start_time.hour,
                                   self.start_time.minute, self.start_time.second)
        td = self.time_step_output
        time_array = [self.start_time + td*i for i in
                      range(self.steps_output)]
        time_array_relative = [td*i for i in range(self.steps_output)]
        return time_array, time_array_relative

    def get_property(self, propname):
        """Get property from history, sorted by status."""
        index_of_first, index_of_last = \
            self.index_of_activation_and_deactivation()
        prop = self.history[propname].copy()
        status = self.history['status'].copy()
        j = np.arange(status.shape[1])
        # Fill arrays with last value before deactivation
        for i in range(status.shape[0]):
            status[i, j > index_of_last[i]] = status[i, index_of_last[i]]
            prop[i, j > index_of_last[i]] = prop[i, index_of_last[i]]

        return prop.T, status.T

    def update_positions(self, x_vel, y_vel):
        """Move particles according to given velocity components.

        This method shall account for projection metrics (a distance
        on a map projection does not necessarily correspond to the same
        distance over true ground (not yet implemented).

        Arguments:
            x_vel and v_vel: floats, velocities in m/s of particle along
                             x- and y-axes of the inherit SRS (proj4).
        """

        geod = pyproj.Geod(ellps='WGS84')

        azimuth = np.degrees(np.arctan2(x_vel, y_vel))  # Direction of motion
        velocity = np.sqrt(x_vel**2 + y_vel**2)  # Velocity in m/s
        velocity = velocity * self.elements.moving  # Do not move frosen elements

        if not self.proj.crs.is_geographic:  # Need to rotate SRS
            # Calculate x,y from lon,lat
            self.elements.x, self.elements.y = self.lonlat2xy(
                self.elements.lon, self.elements.lat)
            # Calculate azimuth orientation of y-axis at particle locations
            delta_y = 1000  # Using delta of 1000 m to calculate azimuth
            lon2, lat2 = self.xy2lonlat(self.elements.x,
                                        self.elements.y + delta_y)
            azimuth_srs = geod.inv(self.elements.lon, self.elements.lat,
                                   lon2, lat2)[0]
            azimuth = azimuth + azimuth_srs

        # Calculate new positions
        self.elements.lon, self.elements.lat, back_az = geod.fwd(
            self.elements.lon, self.elements.lat,
            azimuth, velocity*self.time_step.total_seconds())

        # Check that new positions are valid
        if (self.elements.lon.min() < -180) or (
                self.elements.lon.min() > 360) or (
                self.elements.lat.min() < -90) or (
                self.elements.lat.max() > 90):
            self.logger.info('Invalid new coordinates:')
            self.logger.info(self.elements)
            sys.exit('Quitting')

    def __repr__(self):
        """String representation providing overview of model status."""
        outStr = '===========================\n'
        if hasattr(self, 'history'):
            outStr += self.performance()
            outStr += '===========================\n'
        outStr += 'Model:\t' + type(self).__name__ + \
            '     (OpenDrift version %s)\n' % opendrift.__version__
        outStr += '\t%s active %s particles  (%s deactivated, %s scheduled)\n'\
            % (self.num_elements_active(), self.ElementType.__name__,
               self.num_elements_deactivated(), self.num_elements_scheduled())
        outStr += 'Projection: %s\n' % self.proj4
        variable_groups, reader_groups, missing = self.get_reader_groups()
        outStr += '-------------------\n'
        outStr += 'Environment variables:\n'
        for i, variableGroup in enumerate(variable_groups):
            outStr += '  -----\n'
            readerGroup = reader_groups[i]
            for variable in sorted(variableGroup):
                outStr += '  ' + variable + '\n'
            for i, reader in enumerate(readerGroup):
                outStr += '     ' + str(i+1) + ') ' + reader + '\n'
        if len(self.missing_variables()) > 0:
            outStr += '  -----\n'
            outStr += 'Readers not added for the following variables:\n'
            for variable in sorted(self.missing_variables()):
                outStr += '  ' + variable + '\n'

        lazy_readers = [r for r in self.readers
                            if self.readers[r].is_lazy is True]
        if len(lazy_readers) > 0:
            outStr += '---\nLazy readers:\n'
            for lr in lazy_readers:
                outStr += '  ' + lr + '\n'
        if hasattr(self, 'discarded_readers'):
            outStr += '\nDiscarded readers:\n'
            for dr in self.discarded_readers:
                outStr += '  ' + dr + '\n'
        if hasattr(self, 'time'):
            outStr += '\nTime:\n'
            outStr += '\tStart: %s\n' % (self.start_time)
            outStr += '\tPresent: %s\n' % (self.time)
            if hasattr(self, 'time_step'):
                outStr += '\tCalculation steps: %i * %s - total time: %s\n' % (
                    self.steps_calculation, self.time_step,
                    self.time-self.start_time)
                outStr += '\tOutput steps: %i * %s\n' % (
                    self.steps_output, self.time_step_output)
        if hasattr(self, 'messages'):
            outStr += '-------------------\n'
            outStr += self.get_messages()
        outStr += '===========================\n'
        return outStr

    def store_message(self, message):
        """Store important messages to be displayed to user at end."""
        if not hasattr(self, 'messages'):
            self.messages = []
        self.messages.append(message)

    def get_messages(self):
        """Report any messages stored during simulation."""

        if hasattr(self, 'messages'):
            return str(self.messages).strip('[]') + '\n'
        else:
            return ''

    def add_halo_readers(self):
        """Adding some Thredds and file readers in prioritised order"""

        self.add_readers_from_file(self.test_data_folder() +
            '../../opendrift/scripts/data_sources.txt')

    def calculate_ftle(self, reader=None, delta=None, domain=None,
                       time=None, time_step=None, duration=None, z=0,
                       RLCS=True, ALCS=True):

        if reader is None:
            self.logger.info('No reader provided, using first available:')
            reader = list(self.readers.items())[0][1]
            self.logger.info(reader.name)
        if isinstance(reader, pyproj.Proj):
            proj = reader
        elif isinstance(reader,str):
            proj = pyproj.Proj(reader)
        else:
            proj = reader.proj

        import scipy.ndimage as ndimage
        from opendrift.models.physics_methods import ftle

        if not isinstance(duration, timedelta):
            duration = timedelta(seconds=duration)

        if domain==None:
            xs = np.arange(reader.xmin, reader.xmax, delta)
            ys = np.arange(reader.ymin, reader.ymax, delta)
        else:
            xmin, xmax, ymin, ymax = domain
            xs = np.arange(xmin, xmax, delta)
            ys = np.arange(ymin, ymax, delta)

        X, Y = np.meshgrid(xs, ys)
        lons, lats = proj(X, Y, inverse=True)

        if time is None:
            time = reader.start_time
        if not isinstance(time, list):
            time = [time]
        # dictionary to hold LCS calculation
        lcs = {'time': time, 'lon': lons, 'lat':lats}
        lcs['RLCS'] = np.zeros((len(time), len(ys), len(xs)))
        lcs['ALCS'] = np.zeros((len(time), len(ys), len(xs)))
        T = np.abs(duration.total_seconds())
        for i, t in enumerate(time):
            self.logger.info('Calculating LCS for ' + str(t))
            # Forwards
            if RLCS is True:
                self.reset()
                self.seed_elements(lons.ravel(), lats.ravel(),
                                   time=t, z=z)
                self.run(duration=duration, time_step=time_step)
                f_x1, f_y1 = proj(
                    self.history['lon'].T[-1].reshape(X.shape),
                    self.history['lat'].T[-1].reshape(X.shape))
                lcs['RLCS'][i,:,:] = ftle(f_x1-X, f_y1-Y, delta, T)
            # Backwards
            if ALCS is True:
                self.reset()
                self.seed_elements(lons.ravel(), lats.ravel(),
                                   time=t+duration, z=z)
                self.run(duration=duration, time_step=-time_step)
                b_x1, b_y1 = proj(
                    self.history['lon'].T[-1].reshape(X.shape),
                    self.history['lat'].T[-1].reshape(X.shape))
                lcs['ALCS'][i,:,:] = ftle(b_x1-X, b_y1-Y, delta, T)

        lcs['RLCS'] = np.ma.masked_invalid(lcs['RLCS'])
        lcs['ALCS'] = np.ma.masked_invalid(lcs['ALCS'])
        # Flipping ALCS left-right. Not sure why this is needed
        lcs['ALCS'] = lcs['ALCS'][:,::-1,::-1]

        return lcs

    def center_of_gravity(self, onlysurface=False):
        """
        calculate center of mass and variance of all elements
        returns  (lon,lat), variance
        where (lon,lat) are the coordinates of the center of mass as
        function of time"""
        #lon,lat = self.get_property('lon')[0], self.get_property('lat')[0]
        lon,lat = self.history['lon'], self.history['lat']
        x,y = self.proj(lon, lat)
        if onlysurface==True:
            z = self.history['z']
            submerged = z < 0
            x = np.ma.array(x, mask=submerged)
            y = np.ma.array(y, mask=submerged)
        # center of gravity:
        x_m, y_m = np.ma.mean(x, axis=0), np.ma.mean(y, axis=0)
        center =  self.proj(x_m, y_m, inverse=True)
        one = np.ones_like(x)
        # variance:
        variance = np.ma.mean((x-x_m*one)**2 + (y-y_m*one)**2, axis=0)

        return center,variance

    def reset(self):
        """Preparing OpenDrift object for new run"""
        if not hasattr(self, 'start_time'):
            self.logger.info('Nothing to reset')
            return

        for attr in ['start_time', 'history', 'elements']:
            if hasattr(self, attr):
                delattr(self, attr)
        #del self.start_time
        #del self.history
        #del self.elements
        self.elements_deactivated = self.ElementType()  # Empty array
        self.elements = self.ElementType()  # Empty array
