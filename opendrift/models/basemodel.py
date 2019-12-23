"""
"""
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
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import sys
import os
import glob
import types
import traceback
import logging; logging.captureWarnings(True)
import warnings
from datetime import datetime, timedelta
from collections import OrderedDict
from abc import ABCMeta, abstractmethod, abstractproperty
import netCDF4
from future.utils import iteritems

import numpy as np
import scipy
import pyproj
import configobj, validate
try:
    import matplotlib
    matplotlib.rcParams['legend.numpoints'] = 1
    matplotlib.rcParams['legend.scatterpoints'] = 1
    if ('DISPLAY' not in os.environ and
            'PYCHARM_HOSTED' not in os.environ and
            os.name != 'nt'):
        logging.info('No display found. Using non-interactive Agg backend')
        matplotlib.use('agg')
    import matplotlib.pyplot as plt
    from matplotlib import animation
    from matplotlib.patches import Polygon
    from matplotlib.path import Path
    import cartopy
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except ImportError:
    print('matplotlib and/or cartopy is not available, can not make plots')

import opendrift
from opendrift.readers.basereader import BaseReader, vector_pairs_xy
from opendrift.readers import reader_from_url
from opendrift.models.physics_methods import PhysicsMethods

try:
    basestring
except NameError:
    basestring = str


class OpenDriftSimulation(PhysicsMethods):
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

    configspec_basemodel = '''
            [general]
                use_auto_landmask = boolean(default=True)
                coastline_action = option('none', 'stranding', 'previous', default='stranding')
                time_step_minutes = integer(min=1, max=1440, default=60)
                time_step_output_minutes = integer(min=1, max=1440, default=None)
            [seed]
                ocean_only = boolean(default=True)
            [drift]
                max_age_seconds = float(min=0, default=None)
                scheme = option('euler', 'runge-kutta', 'runge-kutta4', default='euler')
                stokes_drift = boolean(default=True)
                current_uncertainty = float(min=0, max=5, default=0)
                current_uncertainty_uniform = float(min=0, max=5, default=0)
                wind_uncertainty = float(min=0, max=5, default=0)
                relative_wind = boolean(default=False)
                lift_to_seafloor = boolean(default=True)
                truncate_ocean_model_below_m = float(min=0, max=10000, default=None)
                deactivate_north_of = float(min=-90, max=90, default=None)
                deactivate_south_of = float(min=-90, max=90, default=None)
                deactivate_east_of = float(min=-360, max=360, default=None)
                deactivate_west_of = float(min=-360, max=360, default=None)
                use_tabularised_stokes_drift = boolean(default=False)
                tabularised_stokes_drift_fetch = option(5000, 25000, 50000, default=25000)'''

    max_speed = 1  # Assumed max average speed of any element
    required_profiles = None  # Optional possibility to get vertical profiles
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

        # Set default configuration
        self._add_configstring(self.configspec_basemodel)
        self.show_continuous_performance = False

        # Dict to store readers
        self.readers = OrderedDict()  # Dictionary, key=name, value=reader object
        self.priority_list = OrderedDict()
        self.use_block = True  # Set to False if interpolation left to reader

        if not hasattr(self, 'fallback_values'):
            self.fallback_values = {}

        # Make copies of dictionaries so that they are private to each instance
        self.status_categories = ['active']  # Particles are active by default
        self.fallback_values = self.fallback_values.copy()
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
                handler = logging.FileHandler(logfile)
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

        self.timer_start('total time')
        self.timer_start('configuration')

        self.add_metadata('opendrift_version', opendrift.__version__)
        self.logger.info('OpenDriftSimulation initialised (version %s)' %
                     opendrift.__version__)

    def _add_config(self, key, value, comment, overwrite=False):
        """Add a configuration item to the current model."""
        if isinstance(value, list):
            value = 'option(%s, default=\'%s\')' % (
                     str(value)[1:-1], value[0])
        cs = ''  # configstring
        for i, s in enumerate(key.split(':')):
            if i < len(key.split(':')) - 1:
                cs += '\n' + '['*(i+1) + s + ']'*(i+1)
            else:
                cs += '\n\t' + s + ' = ' + value
        newconf = configobj.ConfigObj(configspec=cs.split('\n'),
                    encoding='utf8')
        newconf.validate(validate.Validator())
        if not hasattr(self, 'configobj'):
            self.configobj = newconf
        else:
            self._merge_config(newconf, overwrite=True)

    def _merge_config(self, newconf, overwrite=False):
        if overwrite is True:
            self.configobj.merge(newconf)
            self.configobj.configspec.merge(newconf.configspec)
        else:
            newconf.configspec.merge(self.configobj.configspec)
            self.configobj.configspec = newconf.configspec
            newconf.merge(self.configobj)
            self.configobj = newconf

    def _add_configstring(self, configstring):
        """Add several configuration items from string (INI-format)."""
        newconf = configobj.ConfigObj(configspec=configstring.split('\n'),                      encoding='utf8')
        newconf.validate(validate.Validator())
        if not hasattr(self, 'configobj'):
            self.configobj = newconf
        else:
            self.configobj.merge(newconf)
            self.configobj.configspec.merge(newconf.configspec)

    def set_config(self, key, value):
        """Provide a value for a configuration setting"""
        if key == 'general:use_basemap_landmask':
            warnings.warn("use_basemap_landmask is deprecated in favor of use_auto_landmask", DeprecationWarning)
            key = 'general:use_auto_landmask'

        d = self.configobj
        ds = self.configobj.configspec
        for i, s in enumerate(key.split(':')):
            if s not in d:
                self.list_configspec()
                raise ValueError('Wrong configuration')
            if not isinstance(d[s], dict):
                if ds[s][0:5] == 'float' and value is not None:
                    value = float(value)
                if ds[s][0:5] == 'integ' and value is not None:
                    value = int(value)
                d[s] = value
            else:
                d = d[s]
                ds = ds[s]
        valid = self.configobj.validate(validate.Validator())
        if self.configobj.validate(validate.Validator()) is not True:
            if hasattr(self, 'oiltypes') or hasattr(self, 'leewayprop'):
                # Provide some suggestion if config is very long list (e.g. oiltype/leewaycategory)
                if hasattr(self, 'oiltypes'):
                    elements = self.oiltypes
                elif hasattr(self, 'leewayprop'):
                    elements = [self.leewayprop[la]['OBJKEY'] for la in self.leewayprop]
                import difflib
                matches = difflib.get_close_matches(value, elements, n=20, cutoff=.3)
                containing = [el for el in elements if value in el]
                matches = list(set(matches) | set(containing))
                if len(matches) > 0:
                    matches.sort()
                    suggestion = '\nDid you mean any of these?\n%s' % str(matches)
            else:
                suggestion = ''
            raise ValueError('Wrong configuration:\n\t%s = %s%s' % (s, ds[s], suggestion))

    def _config_hash_to_dict(self, hashstring):
        """Return configobj-dict from section:section:key format"""
        c = self.configobj
        cs = self.configobj.configspec
        for i, s in enumerate(hashstring.split(':')):
            if isinstance(c[s], dict):
                c = c[s]
                cs = cs[s]
        return c, cs, s

    def get_config(self, key):
        """Get value of a configuration setting"""
        c, cs, s = self._config_hash_to_dict(key)
        return c[s]

    def get_configspec(self, key):
        """Get configuration specification of a configuration setting"""
        c, cs, s = self._config_hash_to_dict(key)
        return cs[s]

    def get_configspec_default_min_max(self, key):
        """Get default, min and max config values for configuration setting"""
        c, cs, s = self._config_hash_to_dict(key)
        v = validate.Validator()
        default = v.get_default_value(cs[s])
        datatype = cs[s].split('(')[0]
        if datatype in ['float', 'integer']:
            minimum = cs[s].split('min=')[-1].split(',')[0]
            maximum = cs[s].split('max=')[-1].split(',')[0]
            if datatype == 'float':
                minimum = float(minimum)
                maximum = float(maximum)
            elif datatype == 'integer':
                minimum = int(minimum)
                maximum = int(maximum)
        else:
            minimum = maximum = None
        if datatype == 'option':
            # Parsing options list
            options = cs[s][len('options'):-1]
            options = options.split('default=')[0][:-2]
            import ast
            options = ast.literal_eval('[' + options + ']')
        else:
            options = None
        return default, minimum, maximum, options

    def get_seed_config(self):
        """Get dictionary with all config options under category seed"""
        sc = {}
        for c in self._config_hashstrings():
            if c[0:5] == 'seed:':
                name = c[5:]
                sc[name] = {}
                sn = sc[name]
                sn['default'], sn['min'], sn['max'], sn['options'] = \
                    self.get_configspec_default_min_max(c)
        return sc

    def _config_hashstrings(self):
        """Get list of strings section:section:key for all config items"""
        def fun(k, d, pre):
            path = '%s:%s' % (pre, k) if pre else k
            return path if type(d[k]) is not configobj.Section else \
                ",".join([fun(i,d[k], path) for i in d[k]])
        return ",".join([fun(k,self.configobj, '') for
                         k in self.configobj]).split(',')

    def list_configspec(self):
        """List all possible configuration settings with specifications"""
        keys = self._config_hashstrings()
        str = '\n==============================================\n'
        for key in keys:
            str += '%s [%s] %s\n' % (key, self.get_config(key),
                                   self.get_configspec(key))
        str += '==============================================\n'
        self.logger.info(str)

    def list_config(self):
        """List all possible configuration settings with values"""
        keys = self._config_hashstrings()
        str = '\n=============================================\n'
        for key in keys:
            str += '%s [%s]\n' % (key, self.get_config(key))
        str += '=============================================\n'
        self.logger.info(str)

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
                    (self.elements.ID >= self.newly_seeded_IDs[0]),
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
        warnings.warn ("Use `test_data` fixture", DeprecationWarning)
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
        for category, time in iteritems(self.timing):
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

    def add_reader(self, readers, variables=None):
        """Add one or more readers providing variables used by this model.

        Method may be called subsequently to add more readers
        for other variables.

        Args:
            readers: one or more (list) Reader objects.
            variables: optional, list of strings of standard_name of
                variables to be provided by this/these reader(s).
        """

        # Convert any strings to lists, for looping
        if isinstance(variables, basestring):
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

        if isinstance(urls, basestring):
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
            variables = self.required_variables
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
                    - lists of one ReaderBlock per variable group:
                        - time, x, y, [vars]
        Returns:
            environment: recarray with variables as named attributes,
                         interpolated to requested positions/time.

        '''
        self.timer_start('main loop:readers')
        # Initialise ndarray to hold environment variables
        dtype = [(var, np.float32) for var in variables]
        env = np.ma.array(np.zeros(len(lon))*np.nan, dtype=dtype)

        # Discard any existing readers which are not relevant
        self.discard_irrelevant_readers()

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
            if (self.fallback_values is not None
                    and variable in self.fallback_values):
                env[variable] = np.ma.ones(env[variable].shape)\
                    * self.fallback_values[variable]

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
        if self.get_config('drift:use_tabularised_stokes_drift') is True:
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
        if not isinstance(time, list):
            time = [time]*len(elements)  # Convert to array of same length
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
                loglevel=self.logger.getLogger().getEffectiveLevel())
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
        if tmp_reader is True:
            plt.close()

        return lon, lat

    def seed_elements(self, lon, lat, radius=0, number=None, time=None,
                      cone=False, **kwargs):
        """Seed a given number of particles around given position(s).

        Arguments:
            lon: scalar or array, central longitude(s).
            lat: scalar or array, central latitude(s).
            radius: scalar or array, radius in meters around each lon-lat
                pair, within which particles will be randomly seeded.
            number: integer, total number of particles to be seeded
                Elements are spread equally among the given lon/lat points.
                Default is one particle for each lon-lat pair.
            time: datenum, the time at which particles are seeded/released.
                If time is an array with two elements, elements are seeded
                continously from start/first to end/last time.
            cone: boolean or integer. If True, lon and lat must be two element
                arrays, interpreted as the start and end position of a cone
                within which elements will be seeded. Radius may also be a
                two element array specifying the radius around the points.
            kwargs: keyword arguments containing properties/attributes and
                values corresponding to the actual particle type (ElementType).
                These are forwarded to the ElementType class. All properties
                for which there are no default value must be specified.
        """

        if time is None:
            raise ValueError('Time of seeding must be specified')

        #################################################################
        # Make arrays of all input parameters, with one element per
        # lon/lat pair, for sequential seeding
        #################################################################
        lon = np.atleast_1d(lon).ravel()
        lat = np.atleast_1d(lat).ravel()
        num_points = len(lon)  # Number of lon/lat pairs
        if number is not None and number < num_points:
            raise ValueError('Number of elements must be greater or equal '
                             'to number of points.')

        if num_points == 1 and number is None:
            try:
                number = len(time)
            except:
                number = 1

        if num_points > 1:
            ###############################
            # lon and lat are arrays
            ###############################
            radius_array = np.atleast_1d(radius).ravel()
            if len(radius_array) == 1:
                # If scalar radius is given, apply to all points
                radius_array = radius_array*np.ones(num_points)
            if number is None:
                # Default is one element per given position
                #number = 1*np.ones(num_points)
                number = num_points  # temporarily scalar, should be array
            number = np.atleast_1d(number)
            if len(number) == 1 and num_points > 1:
                number_array = np.floor(number/num_points)*np.ones(num_points)
                # Add remaining elements to first point
                remainder = number - np.floor(number/num_points)*num_points
                number_array[0] = number_array[0] + remainder

            if isinstance(time, datetime):
                time = [time, time]

            if type(time) == list and len(time) == 2:
                td = (time[1]-time[0])/(number-1)  # timestep between points
                if len(td) == 1:
                    td = td[0]
                time_array = [time[0] + i*td for i in range(number[0])]
                indx_time_end = np.cumsum(number_array, dtype=int)
                indx_time_start = np.append([0], indx_time_end[0:-1])
                time_array2 = [time_array[int(indx_time_start[i]):
                                          int(indx_time_end[i])]
                               for i in range(num_points)]
                time_array = time_array2  # Subset of times for this point
            if type(time) == list and len(time) > 2:
                time_array = time

            if cone is True:
                ###################################################
                # lon and lat are start and end points of a cone
                ###################################################
                if len(lon) != 2 or len(lat) != 2:
                    raise ValueError('When cone is True, lon and lat must '
                                     'be 2-element arrays.')

                geod = pyproj.Geod(ellps='WGS84')
                conelonlats = geod.npts(lon[0], lat[0], lon[1], lat[1],
                                        number, radians=False)
                # Seed cone recursively
                lon, lat = list(zip(*conelonlats))
                lon = np.atleast_1d(lon)
                lat = np.atleast_1d(lat)
                number = int(number)
                if len(radius_array) == 1:
                    radius_array = [radius, radius]
                radius_array = np.linspace(radius_array[0], radius_array[1],
                                           number)
                number_array = np.ones(number)
                time_array = [time[0] + i*td for i in np.arange(number)]

            if 'z' in kwargs and isinstance(kwargs['z'], basestring) \
                    and kwargs['z'][0:8] == 'seafloor':
                # We need to fetch seafloor depth from reader
                if ('sea_floor_depth_below_sea_level' not in self.priority_list
                    ) and len(self._lazy_readers()) == 0:
                    raise ValueError('A reader providing the variable '
                                     'sea_floor_depth_below_sea_level must be '
                                     'added before seeding elements at seafloor.')
                if type(time) is list:
                    t = time[0]
                else:
                    t = time
                if not hasattr(self, 'time'):
                    self.time = t
                env, env_profiles, missing = \
                    self.get_environment(['sea_floor_depth_below_sea_level'],
                                         time=t, lon=lon, lat=lat,
                                         z=0*lon, profiles=None)
                # Add M meters if given as 'seafloor+M'
                if len(kwargs['z']) > 8 and kwargs['z'][8] == '+':
                    meters_above_seafloor = np.float(kwargs['z'][9::])
                    self.logger.info('Seeding elements %f meters above seafloor'
                                 % meters_above_seafloor)
                else:
                    meters_above_seafloor = 0
                kwargs['z'] = \
                    -env['sea_floor_depth_below_sea_level'].astype('float32') + meters_above_seafloor

            # Recursively seeding elements around each point
            scalarargs = {}
            for i in range(len(lon)):
                for kwarg in kwargs:
                    try:
                        scalarargs[kwarg] = kwargs[kwarg][i]
                    except:
                        scalarargs[kwarg] = kwargs[kwarg]
                # Make sure to call seed function of base class,
                # not of a specific Children class
                OpenDriftSimulation.seed_elements(
                    self, lon=[lon[i]], lat=[lat[i]], radius=radius_array[i],
                    number=int(number_array[i]),
                    time=time_array[i], cone=False, **scalarargs)
            return

        # Below we have only for single points
        if isinstance(time, datetime) and number > 1:
            time = [time, time]

        if type(time) == list and len(time) == 2 and number > 1:
            td = (time[1]-time[0])/(number-1)  # timestep between points
            # Spread seeding times equally between start/first and end/last
            time_array = [time[0] + i*td for i in range(number)]
        else:
            time_array = time

        geod = pyproj.Geod(ellps='WGS84')
        ones = np.ones(number)
        x = np.random.randn(number)*radius
        y = np.random.randn(number)*radius
        az = np.degrees(np.arctan2(x, y))
        dist = np.sqrt(x*x+y*y)
        kwargs['lon'], kwargs['lat'], az = \
            geod.fwd(lon*ones, lat*ones, az, dist, radians=False)

        if 'z' in kwargs and isinstance(kwargs['z'], basestring) \
                and kwargs['z'][0:8] == 'seafloor':
            # We need to fetch seafloor depth from reader
            # Unfortunately, this is duplication of code above
            if ('sea_floor_depth_below_sea_level' not in self.priority_list
                    ) and len(self._lazy_readers()) == 0:
                raise ValueError('A reader providing the variable '
                                 'sea_floor_depth_below_sea_level must be '
                                 'added before seeding elements at seafloor.')
            if type(time) is list:
                t = time[0]
            else:
                t = time
            if not hasattr(self, 'time'):
                self.time = t
            env, env_profiles, missing = \
                self.get_environment(['sea_floor_depth_below_sea_level'],
                                     t, kwargs['lon'], kwargs['lat'],
                                     0.*ones, None)
            # Add M meters if given as 'seafloor+M'
            if len(kwargs['z']) > 8 and kwargs['z'][8] == '+':
                meters_above_seafloor = np.float(kwargs['z'][9::])
                self.logger.info('Seeding elements %f meters above seafloor'
                             % meters_above_seafloor)
            else:
                meters_above_seafloor = 0
            kwargs['z'] = \
                -env['sea_floor_depth_below_sea_level'].astype('float32') + meters_above_seafloor

        elements = self.ElementType(**kwargs)

        self.schedule_elements(elements, time_array)

    def seed_within_polygon(self, lons, lats, number, **kwargs):
        """Seed a number of elements within given polygon.

        Arguments:
            lon: array of longitudes
            lat: array of latitudes
            number: int, number of elements to be seeded
            kwargs: keyword arguments containing properties/attributes and
                values corresponding to the actual particle type (ElementType).
                These are forwarded to method seed_elements(). All properties
                for which there are no default value must be specified.

        """
        if number == 0:
            return

        lons = np.asarray(lons)
        lats = np.asarray(lats)
        if len(lons) < 3:
            self.logger.info('At least three points needed to make a polygon')
            return
        if len(lons) != len(lats):
            raise ValueError('lon and lat arrays must have same length.')
        poly = Polygon(list(zip(lons, lats)), closed=True)
        # Place N points within the polygons
        proj = pyproj.Proj('+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f '
                           '+lon_0=%f +R=6370997.0 +units=m +ellps=WGS84'
                           % (lats.min(), lats.max(),
                              (lats.min()+lats.max())/2,
                              (lons.min()+lons.max())/2))
        lonlat = poly.get_xy()
        lon = lonlat[:, 0]
        lat = lonlat[:, 1]
        x, y = proj(lon, lat)
        area = 0.0
        for i in range(-1, len(x)-1):
            area += x[i] * (y[i+1] - y[i-1])
        area = abs(area) / 2

        # Make points, evenly distributed
        deltax = np.sqrt(area/number)
        lonpoints = np.array([])
        latpoints = np.array([])
        lonlat = poly.get_xy()
        lon = lonlat[:, 0]
        lat = lonlat[:, 1]
        x, y = proj(lon, lat)
        xvec = np.linspace(x.min() + deltax/2, x.max() - deltax/2,
                           int((x.max()-x.min())/deltax))
        yvec = np.linspace(y.min() + deltax/2, y.max() - deltax/2,
                           int((y.max()-y.min())/deltax))
        x, y = np.meshgrid(xvec, yvec)
        lon, lat = proj(x, y, inverse=True)
        lon = lon.ravel()
        lat = lat.ravel()
        points = np.c_[lon, lat]
        ind = Path(poly.xy).contains_points(points)
        if not any(ind):  # No elements are inside, we seed on border
            lonpoints = np.append(lonpoints, lons[0:number])
            latpoints = np.append(latpoints, lats[0:number])
        else:
            lonpoints = np.append(lonpoints, lon[ind])
            latpoints = np.append(latpoints, lat[ind])
        if len(ind) == 0:
            self.logger.info('Small or irregular polygon, using center point.')
            lonpoints = np.atleast_1d(np.mean(lons))
            latpoints = np.atleast_1d(np.mean(lats))
        # Truncate if too many
        # NB: should also repeat some points, if too few
        lonpoints = lonpoints[0:number]
        latpoints = latpoints[0:number]
        if len(lonpoints) < number:
            # If number of positions is smaller than requested,
            # we duplicate the first ones
            missing = number - len(lonpoints)
            lonpoints = np.append(lonpoints, lonpoints[0:missing])
            latpoints = np.append(latpoints, latpoints[0:missing])

        # Finally seed at calculated positions
        self.seed_elements(lonpoints, latpoints, number=number,
                           **kwargs)

    def seed_from_wkt(self, wkt, number, **kwargs):
        """Seeds elements within (multi)polygons from WKT"""

        try:
            import ogr
            import osr
        except Exception as e:
            self.logger.warning(e)
            raise ValueError('OGR library is needed to parse WKT')

        geom = ogr.CreateGeometryFromWkt(wkt)
        total_area = 0
        for i in range(0, geom.GetGeometryCount()):
            g = geom.GetGeometryRef(i)
            total_area += g.GetArea()

        self.logger.info('Total area of all polygons: %s m2' % total_area)
        num_seeded = 0
        for i in range(0, geom.GetGeometryCount()):
            g = geom.GetGeometryRef(i)
            num_elements = np.int(number*g.GetArea()/total_area)
            if i == geom.GetGeometryCount()-1:
                # For the last feature we seed the remaining number,
                # avoiding difference due to rounding:
                num_elements = number - num_seeded
            self.logger.info('\tSeeding %s elements within polygon number %s' %
                         (num_elements, str(i)))
            try:
                g.Transform(coordTrans)
            except:
                pass
            b = g.GetBoundary()
            if b is not None:
                points = b.GetPoints()
                lons = [p[0] for p in points]
                lats = [p[1] for p in points]
            else:
                # Alternative if OGR is not built with GEOS support
                r = g.GetGeometryRef(0)
                lons = [r.GetX(j) for j in range(r.GetPointCount())]
                lats = [r.GetY(j) for j in range(r.GetPointCount())]

            self.seed_within_polygon(lons, lats, num_elements, **kwargs)
            num_seeded += num_elements


    def seed_from_shapefile(self, shapefile, number,
                            layername=None, featurenum=None, **kwargs):
        """Seeds elements within contours read from a shapefile"""

        try:
            import ogr
            import osr
        except Exception as e:
            self.logger.warning(e)
            raise ValueError('OGR library is needed to read shapefiles.')

        if 'timeformat' in kwargs:
            # Recondstructing time from filename, where 'timeformat'
            # is forwarded to datetime.strptime()
            kwargs['time'] = datetime.strptime(os.path.basename(shapefile),
                                               kwargs['timeformat'])
            del kwargs['timeformat']

        num_seeded_before = self.num_elements_scheduled()

        targetSRS = osr.SpatialReference()
        targetSRS.ImportFromEPSG(4326)
        try:
            s = ogr.Open(shapefile)
        except:
            s = shapefile

        for layer in s:
            if layername is not None and layer.GetName() != layername:
                self.logger.info('Skipping layer: ' + layer.GetName())
                continue
            else:
                self.logger.info('Seeding for layer: %s (%s features)' %
                             (layer.GetDescription(), layer.GetFeatureCount()))

            coordTrans = osr.CoordinateTransformation(layer.GetSpatialRef(),
                                                      targetSRS)

            if featurenum is None:
                featurenum = range(1, layer.GetFeatureCount() + 1)
            else:
                featurenum = np.atleast_1d(featurenum)
            if max(featurenum) > layer.GetFeatureCount():
                raise ValueError('Only %s features in layer.' %
                                 layer.GetFeatureCount())

            # Loop first through all features to determine total area
            layer.ResetReading()
            area_srs = osr.SpatialReference()
            area_srs.ImportFromEPSG(3857)
            areaTransform = osr.CoordinateTransformation(layer.GetSpatialRef(), area_srs)

            areas = np.zeros(len(featurenum))
            for i, f in enumerate(featurenum):
                feature = layer.GetFeature(f - 1)  # Note 1-indexing, not 0
                if feature is not None:
                    gom = feature.GetGeometryRef().Clone()
                    gom.Transform(areaTransform)
                    areas[i] = gom.GetArea()

            total_area = np.sum(areas)
            layer.ResetReading()  # Rewind to first layer
            self.logger.info('Total area of all polygons: %s m2' % total_area)
            # Find number of points per polygon
            numbers = np.round(number*areas/total_area).astype(int)
            numbers[numbers.argmax()] += np.int(number-sum(numbers))

            for i, f in enumerate(featurenum):
                feature = layer.GetFeature(f - 1)
                if feature is None:
                    continue
                num_elements = numbers[i]
                geom = feature.GetGeometryRef()
                self.logger.info('\tSeeding %s elements within polygon number %s' %
                             (num_elements, featurenum[i]))
                try:
                    geom.Transform(coordTrans)
                except:
                    pass
                #b = geom.GetBoundary()
                #if b is not None:
                #    points = b.GetPoints()
                #    lons = [p[0] for p in points]
                #    lats = [p[1] for p in points]
                #else:
                # Alternative if OGR is not built with GEOS support
                r = geom.GetGeometryRef(0)
                lons = [r.GetX(j) for j in range(r.GetPointCount())]
                lats = [r.GetY(j) for j in range(r.GetPointCount())]

                self.seed_within_polygon(lons, lats, num_elements, **kwargs)

    def seed_from_ladim(self, ladimfile, roms):
        """Seed elements from ladim *.rls text file: [time, x, y, z, name]"""

        data = np.loadtxt(ladimfile,
            dtype={'names': ('time', 'x', 'y', 'z'),
                   'formats': ('S20', 'f4', 'f4', 'f4')},
            usecols=(0,1,2,3))

        time = [datetime.strptime(t, "%Y-%m-%dT%H")
                for t in data['time']]
        time = np.array(time)

        lon, lat = roms.xy2lonlat(data['x'], data['y'])
        z = -data['z']

        self.logger.info('Seeding %i elements from %s:' % (len(lon), ladimfile))
        self.logger.info('    Lons: %f to %f' % (lon.min(), lon.max()))
        self.logger.info('    Lats: %f to %f' % (lat.min(), lat.max()))
        self.logger.info('    Depths: %f to %f' % (z.min(), z.max()))
        self.logger.info('    Time: %s to %s' % (time.min(), time.max()))
        elements = self.ElementType(lon=lon, lat=lat, z=-z)

        self.schedule_elements(elements, time)


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
            for varname, profiles in iteritems(self.environment_profiles):
                self.logger.debug('remove items from profile for '+varname)
                if varname is not 'z':
                    self.environment_profiles[varname] = \
                        profiles[:, ~indices]
            self.logger.debug('Removed %i values from environment_profiles.' %
                          (np.sum(indices)))
            #if self.num_elements_active() == 0:
            #    raise ValueError('No more active elements.')  # End simulation

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
        # Check that configuration is proper
        validation = self.configobj.validate(validate.Validator())
        if validation is True:
            self.logger.info('Config validation OK')
        else:
            raise ValueError('Configuration error: ' + str(validation))

        if self.num_elements_scheduled() == 0:
            raise ValueError('Please seed elements before starting a run.')
        self.elements = self.ElementType()

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
                    self.get_environment(self.required_variables,
                                         self.time,
                                         self.elements.lon,
                                         self.elements.lat,
                                         self.elements.z,
                                         self.required_profiles)

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

                if self.num_elements_active() == 0:
                    raise ValueError('No more active elements, quitting.')

                # Store location, in case elements shall be moved back
                self.store_present_positions()

                #####################################################
                self.logger.debug('Calling %s.update()' %
                              type(self).__name__)
                self.timer_start('main loop:updating elements')
                self.update()
                self.timer_end('main loop:updating elements')
                #####################################################

                if self.num_elements_active() == 0:
                    raise ValueError('No active elements, quitting simulation')

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

    def set_up_map(self, corners=None, buffer=.1, delta_lat=None, lscale=None, fast=False, **kwargs):
        """Generate Figure instance on which trajectories are plotted.

           provide corners=[lonmin, lonmax, latmin, latmax] for specific map selection"""

        lons, lats = self.get_lonlats()

        # Initialise map
        if corners==None:
            lonmin = np.nanmin(lons) - buffer*2
            lonmax = np.nanmax(lons) + buffer*2
            latmin = np.nanmin(lats) - buffer
            latmax = np.nanmax(lats) + buffer
        else:
            lonmin = corners[0]
            lonmax = corners[1]
            latmin = corners[2]
            latmax = corners[3]

        if fast:
            self.logger.warning("plotting fast. this will make your plots less accurate.")

            import matplotlib.style as mplstyle
            mplstyle.use(['fast'])

            # use a spherical earth
            axis = 57.29577951308232  # something to do with pi
            globe = ccrs.Globe(ellipse = None, semimajor_axis = axis, semiminor_axis = axis)
            crs = ccrs.Mercator(globe = globe)

            if lscale is None:
                lscale = 'c'
        else:
            crs = ccrs.Mercator()
            if lscale is None:
                lscale = 'auto'

        meanlat = (latmin + latmax)/2
        aspect_ratio = np.float(latmax-latmin) / (np.float(lonmax-lonmin))
        aspect_ratio = aspect_ratio / np.cos(np.radians(meanlat))
        if aspect_ratio > 1:
            fig = plt.figure(figsize=(11./aspect_ratio, 11.))
        else:
            fig = plt.figure(figsize=(11., 11.*aspect_ratio))

        fig.set_tight_layout(True)
        ax = fig.add_subplot(111, projection=crs)  # need '111' for Python 2
        ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())

        def show_landmask(landmask):
            maxn = 512.
            dx = (lonmax - lonmin) / maxn
            dy = (latmax - latmin) / maxn
            dx = max(landmask.dx, dx)
            dy = max(landmask.dy, dy)

            x = np.array([lonmin, lonmax])
            y = np.array([latmin, latmax])
            ndx = np.int32(dx / landmask.dx)  # normalized delta
            ndy = np.int32(dy / landmask.dy)

            xm, ym = landmask.invtransform * (x, y)
            xm = xm.astype(np.int32)
            ym = ym.astype(np.int32)

            img = landmask.mask[ym[0]:ym[1]:ndy, xm[0]:xm[1]:ndx]

            from matplotlib import colors
            cmap = colors.ListedColormap(['white', 'gray'])
            ax.imshow(img, extent=[lonmin, lonmax, latmin, latmax],
                      transform=ccrs.PlateCarree(), cmap=cmap)

        import shapely
        from shapely.geometry import box

        if 'land_binary_mask' in self.priority_list and self.priority_list['land_binary_mask'][0] == 'global_landmask' \
           and not self.readers['global_landmask'].skippoly \
           and self.readers['global_landmask'].mask.extent.contains(box(lonmin, latmin, lonmax, latmax)):

            self.logger.debug("Using existing GSHHS shapes..")
            landmask = self.readers['global_landmask'].mask

            if fast:
                show_landmask(landmask)

            else:
                extent = box(lonmin, latmin, lonmax, latmax)
                extent = shapely.prepared.prep(extent)
                polys = [p for p in landmask.polys.geoms if extent.intersects(p)]

                ax.add_geometries(polys,
                        ccrs.PlateCarree(),
                        facecolor=cfeature.COLORS['land'],
                        edgecolor='black')
        else:
            self.logger.debug ("Adding GSHHS shapes..")

            if fast:
                from opendrift_landmask_data import Landmask
                show_landmask(Landmask(skippoly=True))

            else:
                f = cfeature.GSHHSFeature(scale=lscale, levels=[1],
                        facecolor=cfeature.COLORS['land'])
                ax.add_geometries(
                        f.intersecting_geometries([lonmin, lonmax, latmin, latmax]),
                        ccrs.PlateCarree(),
                        facecolor=cfeature.COLORS['land'],
                        edgecolor='black')


        gl = ax.gridlines(ccrs.PlateCarree(), draw_labels = True)
        gl.xlabels_top = False

        try:
            firstlast = np.ma.notmasked_edges(lons, axis=1)
            index_of_first = firstlast[0][1]
            index_of_last = firstlast[1][1]
        except:
            index_of_last = 0

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except:
            pass

        try:  # Maximise figure window size
            mng.resize(*mng.window.maxsize())
        except:
            pass

        return fig, ax, crs, lons.T, lats.T, index_of_first, index_of_last

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

    def animation(self, buffer=.2, corners=None, filename=None, compare=None,
                  background=None, bgalpha=.5, vmin=None, vmax=None, drifter=None,
                  skip=5, scale=10, color=False, clabel=None,
                  colorbar=True, cmap=None, density=False, show_elements=True,
                  show_trajectories=False,
                  density_pixelsize_m=1000, unitfactor=1, lcs=None,
                  surface_only=False, markersize=20,
                  legend=None, legend_loc='best', fps=10, lscale=None, fast=False):
        """Animate last run."""


        if self.num_elements_total() == 0:
            raise ValueError('Please run simulation before animating')

        start_time = datetime.now()
        if cmap is None:
            cmap = 'jet'
        if isinstance(cmap, basestring):
            cmap = matplotlib.cm.get_cmap(cmap)

        if color is False and background is None and lcs is None and density is False:
            colorbar = False

        markercolor = self.plot_comparison_colors[0]

        if isinstance(density, basestring):
            # Density field is weighted by this variable
            # TODO: not yet implemented!
            density_weight = density
            density = True
        else:
            if density is True:
                density_weight = None

        # x, y are longitude, latitude -> i.e. in a PlateCarree CRS
        gcrs = ccrs.PlateCarree()

        def plot_timestep(i):
            """Sub function needed for matplotlib animation."""
            ax.set_title('%s UTC' % times[i])
            if background is not None:
                map_x, map_y, scalar, u_component, v_component = \
                    self.get_map_background(ax, background,
                                            time=times[i])
                # https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data
                bg.set_array(scalar[:-1,:-1].ravel())
                if type(background) is list:
                    bg_quiv.set_UVC(u_component[::skip, ::skip], v_component[::skip, ::skip])

            if lcs is not None:
                map_x, map_y = map(lcs['lon'], lcs['lat'])
                ax.pcolormesh(
                    map_x, map_y, lcs['ALCS'][i,:,:], alpha=bgalpha,
                    vmin=vmin, vmax=vmax, cmap=cmap)

            if density is True:
                # Update density plot
                pm.set_array(H[i,:,:].ravel())

            # Move points
            if show_elements is True:
                points.set_offsets(np.c_[x[i, range(x.shape[1])],
                                         y[i, range(x.shape[1])]])
                points_deactivated.set_offsets(np.c_[
                    x_deactive[index_of_last_deactivated < i],
                    y_deactive[index_of_last_deactivated < i]])
                if color is not False:  # Update colors
                    points.set_array(colorarray[:, i])
                    if isinstance(color, basestring):
                        points_deactivated.set_array(
                            colorarray_deactivated[
                                index_of_last_deactivated < i])

            if drifter is not None:
                from bisect import bisect_left
                ind = np.max(
                    (0, bisect_left(drifter['time'], times[i]) - 1))
                if i < 1 or i >= len(times)-1 or \
                        drifter['time'][ind] < times[i-1] or \
                        drifter['time'][ind] > times[i+1]:
                    # Do not show when outside time interval
                    drifter_pos.set_offsets([])
                else:
                    drifter_pos.set_offsets(
                        np.c_[drifter['x'][ind], drifter['y'][ind]])

            if show_elements is True:
                if compare is not None:
                    for cd in compare_list:
                        cd['points_other'].set_offsets(np.c_[
                            cd['x_other'][range(cd['x_other'].shape[0]), i],
                            cd['y_other'][range(cd['x_other'].shape[0]), i]])
                        cd['points_other_deactivated'].set_offsets(np.c_[
                            cd['x_other_deactive'][
                                cd['index_of_last_deactivated_other'] < i],
                            cd['y_other_deactive'][
                                cd['index_of_last_deactivated_other'] < i]])
                    return points, cd['points_other']
                else:
                    return points

        # Find map coordinates and plot points with empty data
        fig, ax, crs, x, y, index_of_first, index_of_last = \
            self.set_up_map(buffer=buffer,corners=corners, lscale=lscale, fast=fast)

        if surface_only is True:
            z = self.get_property('z')[0]
            x[z<0] = np.nan
            y[z<0] = np.nan

        if show_trajectories is True:
            ax.plot(x, y, color='gray', alpha=.1, transform = gcrs)

        if color is not False:
            if isinstance(color, basestring):
                colorarray = self.get_property(color)[0].T
                colorarray = colorarray*unitfactor
                colorarray_deactivated = \
                    getattr(self.elements_deactivated, color)
            else:
                colorarray = color
            if vmin is None:
                vmin = colorarray.min()
                vmax = colorarray.max()

        if background is not None:
            map_x, map_y, scalar, u_component, v_component = \
                self.get_map_background(ax, background,
                                        time=self.start_time)
            bg = ax.pcolormesh(map_x, map_y, scalar[:-1,:-1], alpha=bgalpha,
                               antialiased=True, linewidth=0.0, rasterized=True,
                               vmin=vmin, vmax=vmax, cmap=cmap, transform = gcrs)
            if type(background) is list:
                bg_quiv = ax.quiver(map_x[::skip, ::skip],
                                    map_y[::skip, ::skip],
                                    u_component[::skip, ::skip],
                                    v_component[::skip, ::skip], scale=scale, transform = gcrs)

        if lcs is not None:
            if vmin is None:
                vmin = lcs['ALCS'].min()
                vmax = lcs['ALCS'].max()
            lcsh = ax.pcolormesh(lcs['lon'], lcs['lat'], lcs['ALCS'][0,:,:],
                                  vmin=vmin, vmax=vmax, cmap=cmap,
                                  transform = gcrs)

        times = self.get_time_array()[0]
        index_of_last_deactivated = \
            index_of_last[self.elements_deactivated.ID-1]
        if legend is None:
            legend = ['']

        if color is False:
            c = markercolor
        else:
            c = []
        points = ax.scatter([], [], c=c, zorder=10,
                             edgecolor='', cmap=cmap, s=markersize,
                             vmin=vmin, vmax=vmax, label=legend[0], transform = gcrs)
        # Plot deactivated elements, with transparency
        points_deactivated = ax.scatter([], [], color=c, zorder=9,
                                         vmin=vmin, vmax=vmax, s=markersize,
                                         edgecolor='', alpha=.3, transform = gcrs)
        x_deactive, y_deactive = (self.elements_deactivated.lon, self.elements_deactivated.lat)

        if compare is not None:
            compare_list = self._get_comparison_xy_for_plots(compare)

            for cn, cd in enumerate(compare_list):
                if legend != ['']:
                    legstr = legend[cn+1]
                else:
                    legstr = None
                cd['points_other'] = \
                    ax.scatter([], [], color=
                                self.plot_comparison_colors[cn+1],
                                s=markersize,
                                label=legstr, zorder=10, transform = gcrs)
                # Plot deactivated elements, with transparency
                cd['points_other_deactivated'] = \
                    ax.scatter([], [], alpha=.3, zorder=9, color=
                                self.plot_comparison_colors[cn+1],
                                s=markersize, transform = gcrs)

            if legend != ['', '']:
                plt.legend(markerscale=2, loc=legend_loc)

        if density is True:
            cmap.set_under('w')
            H, H_submerged, H_stranded, lon_array, lat_array = \
                self.get_density_array(pixelsize_m=density_pixelsize_m,
                                       weight=density_weight)
            H = H + H_submerged + H_stranded
            H = np.ma.masked_where(H==0, H)
            lat_array, lon_array = np.meshgrid(lat_array, lon_array)
            pm = ax.pcolormesh(lon_array, lat_array, H[0,:,:],
                               vmin=0.1, vmax=vmax, cmap=cmap, transform=gcrs)

        if drifter is not None:
            drifter['x'], drifter['y'] = (drifter['lon'], drifter['lat'])
            #map.plot(drifter['x'], drifter['y'])
            drifter_pos = ax.scatter([], [], color='r',
                                      zorder=15, label='Drifter', transform = gcrs)

        if colorbar is True:
            if color is not False:
                if isinstance(color, basestring) or clabel is not None:
                    if clabel is None:
                        clabel = color
                item = points
            elif density is not False:
                item = pm
                if clabel is None:
                    clabel = 'density'
            elif lcs is not None:
                item = lcsh
                if clabel is None:
                    clabel = 'LCS'
            elif background is not None:
                item = bg
                if clabel is None:
                    clabel = background

            cb = fig.colorbar(item, orientation='horizontal', pad=.05, aspect=30, shrink=.8)
            cb.set_label(clabel)
            cb.set_alpha(1)
            cb.draw_all()

        fig.canvas.draw()
        fig.set_tight_layout(True)
        anim = animation.FuncAnimation(
            plt.gcf(), plot_timestep, blit=False,
            frames=x.shape[0], interval=50)


        if filename is not None:
            self._save_animation(anim, filename, fps)
            self.logger.debug('Time to make animation: %s' %
                          (datetime.now()-start_time))
        else:
            try:
                plt.show()
            except AttributeError:
                pass

    def animation_profile(self, filename=None, compare=None,
                          legend=['', ''], markersize=5, fps=20):
        """Animate vertical profile of the last run."""

        def plot_timestep(i):
            """Sub function needed for matplotlib animation."""
            #plt.gcf().gca().set_title(str(i))
            ax.set_title('%s UTC' % times[i])
            points.set_data(x[range(x.shape[0]), i],
                            z[range(x.shape[0]), i])
            points_deactivated.set_data(
                x_deactive[index_of_last_deactivated < i],
                z_deactive[index_of_last_deactivated < i])

            if compare is not None:
                points_other.set_data(x_other[range(x_other.shape[0]), i],
                                      z_other[range(x_other.shape[0]), i])
                points_other_deactivated.set_data(
                    x_other_deactive[index_of_last_deactivated_other < i],
                    z_other_deactive[index_of_last_deactivated_other < i])
                return points, points_other
            else:
                return points

        # Set up plot
        index_of_first, index_of_last = \
            self.index_of_activation_and_deactivation()
        z = self.get_property('z')[0].T
        x = self.get_property('lon')[0].T
        #seafloor_depth = \
        #    -self.get_property('sea_floor_depth_below_sea_level')[0].T
        plt.close()
        fig = plt.figure(figsize=(10, 6.))  # Suitable aspect ratio
        ax = fig.gca()
        plt.xlabel('Longitude [degrees]')
        plt.ylabel('Depth [m]')
        times = self.get_time_array()[0]
        index_of_last_deactivated = \
            index_of_last[self.elements_deactivated.ID-1]
        points = plt.plot([], [], '.k', label=legend[0],
                          markersize=markersize)[0]
        # Plot deactivated elements, with transparency
        points_deactivated = plt.plot([], [], '.k', alpha=.3)[0]
        x_deactive = self.elements_deactivated.lon
        z_deactive = self.elements_deactivated.z

        if compare is not None:
            if type(compare) is str:
                # Other is given as filename
                other = self.__class__(loglevel=0)
                other.io_import_file(compare)
            else:
                # Other is given as an OpenDrift object
                other = compare
            z_other = other.get_property('z')[0].T
            x_other = self.get_property('lon')[0].T
            points_other = plt.plot(x_other[0, 0], z_other[0, 0], '.r',
                                    label=legend[1],
                                    markersize=markersize)[0]
            # Plot deactivated elements, with transparency
            points_other_deactivated = plt.plot([], [], '.r', alpha=.3)[0]
            x_other_deactive = other.elements_deactivated.lon
            z_other_deactive = other.elements_deactivated.z
            firstlast = np.ma.notmasked_edges(x_other, axis=1)
            index_of_last_other = firstlast[1][1]
            index_of_last_deactivated_other = \
                index_of_last_other[other.elements_deactivated.ID-1]
            xmax = np.maximum(x.max(), x_other.max())
            xmin = np.minimum(x.min(), x_other.min())
            zmax = np.maximum(z.max(), z_other.max())
            zmin = np.minimum(z.min(), z_other.min())
        else:
            xmin = x.min()
            xmax = x.max()
            zmin = z.min()
            zmax = z.max()

        # Set figure limits
        sky = (zmax-zmin)*.1  # Sky height is 10% of water depth
        plt.xlim([xmin, xmax])
        plt.ylim([zmin, sky])
        ax.add_patch(plt.Rectangle((xmin, 0), xmax-xmin, sky,
                     color='lightsteelblue'))
        ax.add_patch(plt.Rectangle((xmin, zmin), xmax-xmin,
                     -zmin, color='cornflowerblue'))

        if legend != ['', '']:
            plt.legend(loc=4)

        anim = animation.FuncAnimation(plt.gcf(), plot_timestep, blit=False,
                                       frames=x.shape[1], interval=150)

        if filename is not None:
            self._save_animation(anim, filename, fps)

        else:
            try:
                plt.show()
            except AttributeError:
                pass

    def _get_comparison_xy_for_plots(self, compare):
        if not type(compare) is list:
            compare = [compare]
        compare_list = [{}]*len(compare)
        for cn, comp in enumerate(compare):
            compare_list[cn] = {}
            cd = compare_list[cn]  # pointer to dict with data
            if type(comp) is str:
                # Other is given as filename
                other = self.__class__(loglevel=0)
                other.io_import_file(comp)
            else:
                # Other is given as an OpenDrift object
                other = comp

            # Find map coordinates of comparison simulations
            cd['x_other'], cd['y_other'] = \
                (other.history['lon'].copy(), other.history['lat'].copy())
            cd['x_other_deactive'], cd['y_other_deactive'] = \
                (other.elements_deactivated.lon.copy(),
                    other.elements_deactivated.lat.copy())
            cd['firstlast'] = np.ma.notmasked_edges(
                cd['x_other'], axis=1)
            cd['index_of_last_other'] = cd['firstlast'][1][1]
            cd['index_of_last_deactivated_other'] = \
                cd['index_of_last_other'][other.elements_deactivated.ID-1]

        return compare_list

    def plot(self, background=None, buffer=.2, corners=None, linecolor=None, filename=None,
             show=True, vmin=None, vmax=None, compare=None, cmap='jet',
             lvmin=None, lvmax=None, skip=2, scale=10, show_scalar=True,
             contourlines=False, trajectory_dict=None, colorbar=True,
             linewidth=1, lcs=None, show_particles=True,
             density_pixelsize_m=1000,
             surface_color=None, submerged_color=None, markersize=20,
             title='auto', legend=True, legend_loc='best', lscale=None,
             fast=False, **kwargs):
        """Basic built-in plotting function intended for developing/debugging.

        Plots trajectories of all particles.
        Positions marked with colored stars:
        - green: all start positions
        - red: deactivated particles
        - blue: particles still active at end of simulation

        Requires availability of Cartopy.

        Arguments:
            background: string, name of variable (standard_name) which will
                be plotted as background of trajectories, provided that it
                can be read with one of the available readers.
            buffer: float; spatial buffer of plot in degrees of
                longitude/latitude around particle collection.
            background: name of variable to be plotted as background field.
                        Use two element list for vector fields, e.g.
                        ['x_wind', 'y_wind']
            vmin, vmax: minimum and maximum values for colors of background.
            linecolor: name of variable to be used for coloring trajectories.
            lvmin, lvmax: minimum and maximum values for colors of trajectories.
            lscale (string): resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto'). default is 'auto'.
            fast (bool): use some optimizations to speed up plotting at the cost of accuracy
        """


        mappable = None

        if self.num_elements_total() == 0:
            raise ValueError('Please run simulation before plotting')

        start_time = datetime.now()

        # x, y are longitude, latitude -> i.e. in a PlateCarree CRS
        gcrs = ccrs.PlateCarree()
        fig, ax, crs, x, y, index_of_first, index_of_last = \
            self.set_up_map(buffer=buffer,corners=corners, lscale=lscale, fast=fast, **kwargs)

        markercolor = self.plot_comparison_colors[0]

        # The more elements, the more transparent we make the lines
        min_alpha = 0.1
        max_elements = 5000.0
        alpha = min_alpha**(2*(self.num_elements_total()-1)/(max_elements-1))
        alpha = np.max((min_alpha, alpha))
        if legend is False:
            legend = None
        if hasattr(self, 'history') and linewidth != 0:
            # Plot trajectories
            if linecolor is None:
                if compare is not None and legend is not None:
                    if legend is True:
                        if hasattr(compare, 'len'):
                            numleg = len(compare)
                        else:
                            numleg = 2
                        legend = ['Simulation %d' % (i+1) for i in
                                  range(numleg)]
                    ax.plot(x[:,0], y[:,0], color='gray', alpha=alpha, label=legend[0], linewidth=linewidth, transform = gcrs)
                    ax.plot(x, y, color='gray', alpha=alpha, label='_nolegend_', linewidth=linewidth, transform = gcrs)
                else:
                    ax.plot(x, y, color='gray', alpha=alpha, linewidth=linewidth, transform = gcrs)
            else:
                colorbar = True
                # Color lines according to given parameter
                try:
                    if isinstance(linecolor, basestring):
                        param = self.history[linecolor]
                    else:
                        param = linecolor
                except:
                    raise ValueError(
                        'Available parameters to be used for linecolors: ' +
                        str(self.history.dtype.fields))
                from matplotlib.collections import LineCollection
                for i in range(x.shape[1]):
                    vind = np.arange(index_of_first[i], index_of_last[i] + 1)
                    points = np.array(
                        [x[vind, i].T, y[vind, i].T]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]],
                                              axis=1)
                    if lvmin is None:
                        lvmin = param.min()
                        lvmax = param.max()
                    lc = LineCollection(segments,
                                        #cmap=plt.get_cmap('Spectral'),
                                        cmap=cmap,
                                        norm=plt.Normalize(lvmin, lvmax), transform = gcrs)
                    #lc.set_linewidth(3)
                    lc.set_array(param.T[vind, i])
                    mappable = ax.add_collection(lc)
                #axcb = fig.colorbar(lc, ax = ax, orientation = 'horizontal')
                #try:  # Add unit to colorbar if available
                #    colorbarstring = linecolor + '  [%s]' % \
                #        (self.history_metadata[linecolor]['units'])
                #except:
                #    colorbarstring = linecolor
                ##axcb.set_label(colorbarstring)
                #axcb.set_label(colorbarstring, size=14)
                #axcb.ax.tick_params(labelsize=14)

        if compare is None:
            label_initial = 'initial (%i)' % x.shape[1]
            label_active = 'active (%i)' % (x.shape[1] - self.num_elements_deactivated())
            color_initial = self.status_colors['initial']
            color_active = self.status_colors['active']
        else:
            label_initial = None
            label_active = None
            color_initial = 'gray'
            color_active = 'gray'
        if show_particles is True:
            ax.scatter(x[index_of_first, range(x.shape[1])],
                       y[index_of_first, range(x.shape[1])],
                       s=markersize,
                       zorder=10, edgecolor=markercolor, linewidths=.2,
                       color=color_initial, label=label_initial,
                       transform = gcrs)
            if surface_color is not None:
                color_active = surface_color
                label_active = 'surface'
            ax.scatter(x[index_of_last, range(x.shape[1])],
                       y[index_of_last, range(x.shape[1])],
                       s=markersize, zorder=3,
                       edgecolor=markercolor, linewidths=.2,
                       color=color_active, label=label_active,
                       transform = gcrs)
            #if submerged_color is not None:
            #    map.scatter(x[range(x.shape[0]), index_of_last],
            #                y[range(x.shape[0]), index_of_last], s=markersize,
            #                zorder=3, edgecolor=markercolor, linewidths=.2,
            #                color=submerged_color, label='submerged')

            x_deactivated, y_deactivated = (self.elements_deactivated.lon,
                                               self.elements_deactivated.lat)
            # Plot deactivated elements, labeled by deactivation reason
            for statusnum, status in enumerate(self.status_categories):
                if status == 'active':
                    continue  # plotted above
                if status not in self.status_colors:
                    # If no color specified, pick an unused one
                    for color in ['red', 'blue', 'green', 'black', 'gray',
                                  'cyan', 'DarkSeaGreen', 'brown']:
                        if color not in self.status_colors.values():
                            self.status_colors[status] = color
                            break
                indices = np.where(self.elements_deactivated.status == statusnum)
                if len(indices[0]) > 0:
                    if (status == 'seeded_on_land' or
                        status == 'seeded_at_nodata_position'):
                        zorder = 11
                    else:
                        zorder = 3
                    if compare is not None:
                        legstr = None
                    else:
                        legstr = '%s (%i)' % (status, len(indices[0]))
                    if compare is None:
                        color_status = self.status_colors[status]
                    else:
                        color_status = 'gray'
                    ax.scatter(x_deactivated[indices], y_deactivated[indices],
                                s=markersize,
                                zorder=zorder, edgecolor=markercolor, linewidths=.1,
                                color=color_status, label=legstr,
                                transform = gcrs)

        if compare is not None:
            cd = self._get_comparison_xy_for_plots(compare)
            for i, c in enumerate(cd):
                if legend != None:
                    legstr = legend[i+1]
                else:
                    legstr = None
                ax.plot(c['x_other'].T[:,0], c['y_other'].T[:,0], self.plot_comparison_colors[i+1] + '-', label=legstr, linewidth=linewidth, transform = gcrs)
                ax.plot(c['x_other'].T, c['y_other'].T, self.plot_comparison_colors[i+1] + '-', label='_nolegend_', linewidth=linewidth, transform = gcrs)
                ax.scatter(c['x_other'][range(c['x_other'].shape[0]), c['index_of_last_other']],
                    c['y_other'][range(c['y_other'].shape[0]), c['index_of_last_other']],
                    s=markersize,
                    zorder=3, edgecolor=markercolor, linewidths=.2,
                    color=self.plot_comparison_colors[i+1], transform = gcrs)

        try:
            if legend is not None:# and compare is None:
                plt.legend(loc=legend_loc, markerscale=2)
        except Exception as e:
            self.logger.warning('Cannot plot legend, due to bug in matplotlib:')
            self.logger.warning(traceback.format_exc())

        if background is not None:
            if hasattr(self, 'time'):
                time = self.time - self.time_step_output
            else:
                time = None
            if background == 'residence':
                scalar,lon_res,lat_res = self.get_residence_time(
                    pixelsize_m=density_pixelsize_m)
                scalar[scalar==0] = np.nan
                lon_res, lat_res = np.meshgrid(lon_res[0:-1], lat_res[0:-1])
                lon_res = lon_res.T
                lat_res = lat_res.T
                map_x, map_y = (lon_res, lat_res)
            else:
                map_x, map_y, scalar, u_component, v_component = \
                    self.get_map_background(ax, background, time=time)
                                        #self.time_step_output)

            if show_scalar is True:
                if contourlines is False:
                    scalar = np.ma.masked_invalid(scalar)
                    mappable = ax.pcolormesh(map_x, map_y, scalar, alpha=1,
                                   vmin=vmin, vmax=vmax, cmap=cmap, transform = gcrs)
                else:
                    if contourlines is True:
                        CS = ax.contour(map_x, map_y, scalar,
                                         colors='gray', transform = gcrs)
                    else:
                        # contourlines is an array of values
                        CS = ax.contour(map_x, map_y, scalar, contourlines,
                                         colors='gray', transform = gcrs)
                    plt.clabel(CS, fmt='%g')

        if mappable is not None:
            cb = fig.colorbar(mappable, orientation='horizontal', pad=.05, aspect=30, shrink=.8)
            cb.set_label(str(background))

        if type(background) is list:
            delta_x = (map_x[1,2] - map_x[1,1])/2.
            delta_y = (map_y[2,1] - map_y[1,1])/2.
            ax.quiver(map_x[::skip, ::skip] + delta_x, map_y[::skip, ::skip] + delta_y,
                      u_component[::skip, ::skip],
                      v_component[::skip, ::skip], scale=scale, transform = gcrs)

        if lcs is not None:
            map_x_lcs, map_y_lcs = (lcs['lon'], lcs['lat'])
            ax.pcolormesh(
                map_x_lcs, map_y_lcs, lcs['ALCS'][0,:,:], alpha=1,
                vmin=vmin, vmax=vmax, cmap=cmap, transform = gcrs)

        if title is not None:
            if title is 'auto':
                if hasattr(self, 'time'):
                    plt.title(type(self).__name__ + '  %s to %s UTC (%i steps)' %
                              (self.start_time.strftime('%Y-%m-%d %H:%M'),
                               self.time.strftime('%Y-%m-%d %H:%M'),
                               self.steps_output))
                else:
                    plt.title(type(self).__name__ + ' - %i elements seeded at %s UTC' %
                              (self.num_elements_scheduled(),
                               self.elements_scheduled_time[0].strftime(
                               '%Y-%m-%d %H:%M')))
            else:
                plt.title(title)

        if trajectory_dict is not None:
            self._plot_trajectory_dict(ax, trajectory_dict)

        #plt.gca().tick_params(labelsize=14)

        fig.canvas.draw()
        fig.set_tight_layout(True)
        if filename is not None:
            plt.savefig(filename)
            self.logger.info('Time to make plot: ' +
                         str(datetime.now() - start_time))
        else:
            if show is True:
                plt.show()

        return ax, plt

    def _plot_trajectory_dict(self, ax, trajectory_dict):
        '''Plot provided trajectory along with simulated'''
        time = trajectory_dict['time']
        time = np.array(time)
        i = np.where((time>=self.start_time) & (time<=self.time))[0]
        x, y = (trajectory_dict['lon'][i], trajectory_dict['lat'][i])
        ls = trajectory_dict['linestyle']

        gcrs = ccrs.PlateCarree()

        ax.plot(x, y, ls, linewidth=2, transform = gcrs)
        ax.plot(x[0], y[0], 'ok', transform = gcrs)
        ax.plot(x[-1], y[-1], 'xk', transform = gcrs)

    def get_map_background(self, ax, background, time=None):
        # Get background field for plotting on map or animation
        if type(background) is list:
            variable = background[0]  # A vector is requested
        else:
            variable = background  # A scalar is requested
        for readerName in self.readers:
            reader = self.readers[readerName]
            if variable in reader.variables:
                if time is None or (time>= reader.start_time
                        and time <= reader.end_time) or (
                        reader.always_valid is True):
                    break
        if time is None:
            if hasattr(self, 'elements_scheduled_time'):
                # Using time of first seeded element
                time = self.elements_scheduled_time[0]

        # Get reader coordinates covering given map area
        axisproj = pyproj.Proj(ax.projection.proj4_params)
        xmin, xmax, ymin, ymax = ax.get_extent(ccrs.PlateCarree())
        cornerlons = np.array([xmin, xmin, xmax, xmax])
        cornerlats = np.array([ymin, ymax, ymin, ymax])
        reader_x, reader_y = reader.lonlat2xy(cornerlons, cornerlats)
        reader_x = np.linspace(reader_x.min(), reader_x.max(), 10)
        reader_y = np.linspace(reader_y.min(), reader_y.max(), 10)

        data = reader.get_variables(
            background, time, reader_x, reader_y, None, block=True)
        reader_x, reader_y = np.meshgrid(data['x'], data['y'])
        if type(background) is list:
            u_component = data[background[0]]
            v_component = data[background[1]]
            scalar = np.sqrt(u_component**2 + v_component**2)
            # NB: rotation not completed!
            u_component, v_component = reader.rotate_vectors(
                reader_x, reader_y, u_component, v_component,
                reader.proj, ccrs.PlateCarree(
                globe=ccrs.Globe(datum='WGS84', ellipse='WGS84')
                ).proj4_init)
        else:
            scalar = data[background]
            u_component = v_component = None

        # Shift one pixel for correct plotting
        reader_x = reader_x - reader.delta_x
        reader_y = reader_y - reader.delta_y

        rlons, rlats = reader.xy2lonlat(reader_x, reader_y)
        if rlons.max() > 360:
            rlons = rlons - 360
        map_x, map_y = (rlons, rlats)

        scalar = np.ma.masked_invalid(scalar)
        if hasattr(reader, 'convolve'):
            from scipy import ndimage
            N = reader.convolve
            if isinstance(N, (int, np.integer)):
                kernel = np.ones((N, N))
                kernel = kernel/kernel.sum()
            else:
                kernel = N
            self.logger.debug('Convolving variables with kernel: %s' % kernel)
            scalar = ndimage.convolve(
                scalar, kernel, mode='nearest')

        return map_x, map_y, scalar, u_component, v_component

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

        from netCDF4 import Dataset, date2num
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

    def write_geotiff(self, filename, pixelsize_km=.2):
        '''Write one GeoTiff image per timestep.

        filename should contain date identifiers, e.g. 'img_%Y%m%d_%H%M.tif'
        https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
        '''

        try:
            import gdal
            import osr
        except:
            raise ValueError('GDAL is needed to write geotiff images.')
        import matplotlib.pyplot as plt
        driver = gdal.GetDriverByName('GTiff')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        colortable = gdal.ColorTable()
        colortable.SetColorEntry(0,(255,255,255, 0))
        colortable.SetColorEntry(1,(0,0,0, 255))
        colortable.SetColorEntry(2,(255,0,0, 255))
        colortable.SetColorEntry(3,(0,255,0, 255))
        colortable.SetColorEntry(4,(0,0,255, 255))

        lon = self.get_property('lon')[0]
        lat = self.get_property('lat')[0]
        status = self.get_property('status')[0]
        times = self.get_time_array()[0]
        deltalat = pixelsize_km/111.0  # km to degrees
        deltalon = deltalat/np.cos(np.radians((lat.min() + lat.max())/2))
        lat_array = np.arange(lat.min()-deltalat, lat.max()+deltalat, deltalat)
        lon_array = np.arange(lon.min()-deltalat, lon.max()+deltalon, deltalon)
        ilon = (np.round((lon-lon.min())/deltalon)).astype(int)
        ilat = (np.round((lat-lat.min())/deltalat)).astype(int)
        # Setting masked values to zero, for use as indices
        ilon[ilon.mask] = 0
        ilat[ilat.mask] = 0
        status[ilon.mask] = 0
        image = np.zeros((len(times), len(lon_array),
                          len(lat_array))).astype(int)
        geotransform = [lon_array.min(), deltalon, 0,
                        lat_array.max(), 0, -deltalat]
        for i, t in enumerate(times):
            image[i, ilon[i,:], ilat[i,:]] = status[i, :] + 1
            filename_i = t.strftime(filename)
            ds = driver.Create(filename_i, len(lon_array), len(lat_array),
                               1, gdal.GDT_Byte, )
            ds.SetProjection(srs.ExportToWkt())
            ds.SetGeoTransform(geotransform)
            outband=ds.GetRasterBand(1)
            outband.SetNoDataValue(0)
            outband.WriteArray(np.fliplr(image[i, :, :]).transpose())
            outband.SetColorTable(colortable)
            ds = None

    def get_time_array(self):
        """Return a list of output times of last run."""
        td = self.time_step_output
        time_array = [self.start_time + td*i for i in
                      range(self.steps_output)]
        time_array_relative = [td*i for i in range(self.steps_output)]
        return time_array, time_array_relative

    def plot_environment(self, filename=None):
        """Plot mean wind and current velocities of element of last run."""
        x_wind = self.get_property('x_wind')[0]
        y_wind = self.get_property('y_wind')[0]
        wind = np.sqrt(x_wind**2 + y_wind**2)
        x_sea_water_velocity = self.get_property('x_sea_water_velocity')[0]
        y_sea_water_velocity = self.get_property('y_sea_water_velocity')[0]
        current = np.sqrt(x_sea_water_velocity**2 + y_sea_water_velocity**2)
        wind = np.ma.mean(wind, axis=1)
        current = np.ma.mean(current, axis=1)
        time, time_relative = self.get_time_array()
        time = np.array([t.total_seconds()/3600. for t in time_relative])

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(time, wind, 'b', label='wind speed')
        ax1.set_ylabel('Wind speed  [m/s]', color='b')
        ax1.set_xlim([0, time[-1]])
        ax1.set_ylim([0, wind.max()])

        ax2 = ax1.twinx()
        ax2.plot(time, current, 'r', label='current speed')
        ax2.set_ylabel('Current speed  [m/s]', color='r')
        ax2.set_xlim([0, time[-1]])
        ax2.set_ylim([0, current.max()])
        for tl in ax1.get_yticklabels():
            tl.set_color('b')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        ax1.set_xlabel('Time  [hours]')
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)

    def plot_property(self, prop, mean=False):
        """Basic function to plot time series of any element properties."""
        import matplotlib.pyplot as plt
        from matplotlib import dates

        hfmt = dates.DateFormatter('%d %b %Y %H:%M')
        plt.close()
        fig = plt.figure()
        ax = fig.gca()
        ax.xaxis.set_major_formatter(hfmt)
        plt.xticks(rotation='vertical')
        times = [self.start_time + n*self.time_step_output
                 for n in range(self.steps_output)]
        data = self.history[prop].T[0:len(times), :]
        if mean is True:  # Taking average over elements
            data = np.mean(data, axis=1)
        plt.plot(times, data)
        plt.title(prop)
        plt.xlabel('Time  [UTC]')
        try:
            plt.ylabel('%s  [%s]' %
                       (prop, self.elements.variables[prop]['units']))
        except:
            plt.ylabel(prop)
        plt.subplots_adjust(bottom=.3)
        plt.grid()
        plt.show()

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

    def _save_animation(self, anim, filename, fps):
        self.logger.info('Saving animation to ' + filename + '...')

        import sys
        if 'sphinx_gallery' in sys.modules:
            import os
            adir = os.path.realpath('../docs/source/gallery/animations')
            os.makedirs(adir, exist_ok=True)
            filename = os.path.join(adir, filename)

        try:
            if filename[-4:] == '.gif':  # GIF
                self.logger.info('Making animated gif...')
                try:
                    anim.save(filename, fps=fps, writer='imagemagick')
                except:  # For old version of matplotlib
                    anim.save(filename, fps=fps)#, clear_temp=False)
                    os.system('convert -delay %i _tmp*.png %s' %
                            (np.abs(self.time_step_output.total_seconds())/
                             3600.*24., filename))
                    self.logger.info('Deleting temporary figures...')
                    tmp = glob.glob('_tmp*.png')
                    for tfile in tmp:
                        os.remove(tfile)
            else:  # MP4
                try:
                    try:
                        # For perfect quality, but larger file size
                        #FFwriter=animation.FFMpegWriter(fps=fps, extra_args=['-vcodec', 'libx264'])
                        FFwriter=animation.FFMpegWriter(fps=fps,
                            codec='libx264', bitrate=1800,
                            extra_args=['-profile:v', 'baseline',
                                        '-pix_fmt', 'yuv420p', '-an'])
                        anim.save(filename, writer=FFwriter)
                    except Exception as e:
                        self.logger.info(e)
                        anim.save(filename, fps=fps, bitrate=1800,
                                  extra_args=['-an',
                                              '-pix_fmt', 'yuv420p'])
                except Exception as e:
                    self.logger.info(e)
                    anim.save(filename, fps=fps)
                    self.logger.warning('Animation might not be HTML5 compatible.')

        except Exception as e:
            self.logger.info('Could not save animation:')
            self.logger.info(e)
            self.logger.debug(traceback.format_exc())

        if 'sphinx_gallery' in sys.modules:
            plt.close('all')

    def calculate_ftle(self, reader=None, delta=None,
                       time=None, time_step=None, duration=None,
                       RLCS=True, ALCS=True):

        if reader is None:
            self.logger.info('No reader provided, using first available:')
            reader = list(self.readers.items())[0][1]
            self.logger.info(reader.name)
        import scipy.ndimage as ndimage
        from opendrift.models.physics_methods import ftle

        if not isinstance(duration, timedelta):
            duration = timedelta(seconds=duration)

        xs = np.arange(reader.xmin, reader.xmax, delta)
        ys = np.arange(reader.ymin, reader.ymax, delta)
        X, Y = np.meshgrid(xs, ys)
        lons, lats = reader.xy2lonlat(X, Y)

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
                                   time=t)
                self.run(duration=duration, time_step=time_step)
                f_x1, f_y1 = reader.lonlat2xy(
                    self.history['lon'].T[-1].reshape(X.shape),
                    self.history['lat'].T[-1].reshape(X.shape))
                lcs['RLCS'][i,:,:] = ftle(f_x1-X, f_y1-Y, delta, T)
            # Backwards
            if ALCS is True:
                self.reset()
                self.seed_elements(lons.ravel(), lats.ravel(),
                                   time=t+duration)
                self.run(duration=duration, time_step=-time_step)
                b_x1, b_y1 = reader.lonlat2xy(
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

