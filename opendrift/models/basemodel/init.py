from abc import abstractproperty
import copy
import logging
import numpy as np
import pyproj

from .state import State

logger = logging.getLogger(__name__)

CONFIG_LEVEL_ESSENTIAL = 1
CONFIG_LEVEL_BASIC = 2
CONFIG_LEVEL_ADVANCED = 3


class Init(State):
    _config = None
    origin_marker = None

    def __new__(cls):
        """
        We override __new__ so that this initialization is always done before sub-class __init__.
        """
        i = super().__new__(cls)
        i._config = dict()
        i._add_config({
            # type, default, min, max, enum, important, value, units, description
            'general:use_auto_landmask': {
                'type':
                'bool',
                'default':
                True,
                'description':
                'A built-in GSHHG global landmask is used if True, '
                'otherwise landmask is taken from reader or fallback value.',
                'level':
                CONFIG_LEVEL_ADVANCED
            },
            'general:coastline_action': {
                'type':
                'enum',
                'enum': ['none', 'stranding', 'previous'],
                'default':
                'stranding',
                'level':
                CONFIG_LEVEL_BASIC,
                'description':
                'None means that objects may also move over land. '
                'stranding means that objects are deactivated if they hit land. '
                'previous means that objects will move back to the previous location '
                'if they hit land'
            },
            'general:time_step_minutes': {
                'type':
                'float',
                'min':
                .01,
                'max':
                1440,
                'default':
                60,
                'units':
                'minutes',
                'level':
                CONFIG_LEVEL_BASIC,
                'description':
                'Calculation time step used for the simulation. The output time step may '
                'be equal or larger than this.'
            },
            'general:time_step_output_minutes': {
                'type':
                'float',
                'min':
                1,
                'max':
                1440,
                'default':
                None,
                'units':
                'minutes',
                'level':
                CONFIG_LEVEL_BASIC,
                'description':
                'Output time step, i.e. the interval at which output is saved. This must be larger than '
                'the calculation time step, and be an integer multiple of this.'
            },
            'seed:ocean_only': {
                'type':
                'bool',
                'default':
                True,
                'description':
                'If True, elements seeded on land will be moved to the closest '
                'position in ocean',
                'level':
                CONFIG_LEVEL_ADVANCED
            },
            'seed:number': {
                'type': 'int',
                'default': 1,
                'min': 1,
                'max': 100000000,
                'units': 1,
                'description': 'The number of elements for the simulation.',
                'level': CONFIG_LEVEL_BASIC
            },
            'drift:max_age_seconds': {
                'type': 'float',
                'default': None,
                'min': 0,
                'max': np.inf,
                'units': 'seconds',
                'description':
                'Elements will be deactivated when this age is reached',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:advection_scheme': {
                'type':
                'enum',
                'enum': ['euler', 'runge-kutta', 'runge-kutta4'],
                'default':
                'euler',
                'level':
                CONFIG_LEVEL_ADVANCED,
                'description':
                'Numerical advection scheme for ocean current advection'
            },
            'drift:current_uncertainty': {
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 5,
                'units': 'm/s',
                'description':
                'Add gaussian perturbation with this standard deviation to current components at each time step',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:current_uncertainty_uniform': {
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 5,
                'units': 'm/s',
                'description':
                'Add gaussian perturbation with this standard deviation to current components at each time step',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:horizontal_diffusivity': {
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 100000,
                'units': 'm2/s',
                'description': 'Add horizontal diffusivity (random walk)',
                'level': CONFIG_LEVEL_BASIC
            },
            'drift:wind_uncertainty': {
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 5,
                'units': 'm/s',
                'description':
                'Add gaussian perturbation with this standard deviation to wind components at each time step.',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:relative_wind': {
                'type': 'bool',
                'default': False,
                'description':
                'If True, wind drift is calculated for absolute wind (wind vector minus ocean surface current vector).',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:deactivate_north_of': {
                'type': 'float',
                'default': None,
                'min': -90,
                'max': 90,
                'units': 'degrees',
                'description':
                'Elements are deactivated if the move further north than this limit',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:deactivate_south_of': {
                'type': 'float',
                'default': None,
                'min': -90,
                'max': 90,
                'units': 'degrees',
                'description':
                'Elements are deactivated if the move further south than this limit',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:deactivate_east_of': {
                'type': 'float',
                'default': None,
                'min': -360,
                'max': 360,
                'units': 'degrees',
                'description':
                'Elements are deactivated if the move further east than this limit',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:deactivate_west_of': {
                'type': 'float',
                'default': None,
                'min': -360,
                'max': 360,
                'units': 'degrees',
                'description':
                'Elements are deactivated if the move further west than this limit',
                'level': CONFIG_LEVEL_ADVANCED
            },
        })

        return i

    def get_config(self, key):
        if not key in self._config:
            raise ValueError('No config setting named %s' % key)
        return (self._config[key]['value'])

    def list_config(self, prefix=''):
        """List all possible configuration settings with values"""
        str = '\n=============================================\n'
        for key in self._config:
            if key.startswith(prefix):
                str += '%s [%s]\n' % (key, self.get_config(key))
        str += '=============================================\n'
        logger.info(str)

    def list_configspec(self, prefix=''):
        """Readable formatting of config specification"""
        for c, i in self._config.items():
            if c.startswith(prefix):
                val = i['value'] if 'value' in i else None

                val = val(self) if callable(val) else val

                if i['type'] == 'bool':
                    rang = ''
                elif i['type'] in ['float', 'int']:
                    rang = 'min: %s, max: %s [%s]' % (i['min'], i['max'],
                                                      i['units'])
                elif i['type'] == 'enum':
                    rang = i['enum']
                print('%-35s [%s] %-5s %s %s...' %
                      (c, val, i['type'], rang, i['description'][0:20]))

    def get_configspec(self, prefix='', level=[1, 2, 3]):
        if not isinstance(level, list):
            level = [level]
        configspec = {
            k: v
            for (k, v) in self._config.items()
            if k.startswith(prefix) and self._config[k]['level'] in level
        }
        return configspec

    def set_config(self, key, value):
        if not key in self._config:
            self.list_config()
            raise ValueError('No config setting named %s' % key)
        i = self._config[key]
        if i['type'] == 'bool':
            if value not in [True, False]:
                raise ValueError('Config value %s must be True or False' % key)
        elif i['type'] in ['float', 'int'] and value is not None:
            if (i['min'] is not None
                    and value < i['min']) or (i['max'] is not None
                                              and value > i['max']):
                raise ValueError('Config value %s must be between %s and %s' %
                                 (key, i['min'], i['max']))
            if i['type'] == 'float' and value is not None:
                value = float(value)
            elif i['type'] == 'int' and value is not None:
                value = int(value)
        elif i['type'] == 'enum':
            if value not in i['enum']:
                suggestion = ''
                if len(i['enum']) > 5:
                    import difflib
                    lowercase_list = [s.lower() for s in i['enum']]
                    lowercase_mapping = {
                        lc: oc
                        for lc, oc in zip(lowercase_list, i['enum'])
                    }
                    matches = difflib.get_close_matches(value.lower(),
                                                        lowercase_list,
                                                        n=20,
                                                        cutoff=.3)
                    containing = [
                        e for e in lowercase_list if value.lower() in e
                    ]
                    matches = list(set(matches) | set(containing))
                    if len(matches) > 0:
                        matches = [
                            lowercase_mapping[match] for match in matches
                        ]
                        matches.sort()
                        suggestion = '\nDid you mean any of these?\n%s' % str(
                            matches)
                raise ValueError(
                    'Wrong configuration (%s=%s), possible values are:\n\t%s\n%s'
                    % (key, value, i['enum'], suggestion))

        self._config[key]['value'] = value

    def _set_config_default(self, key, value):
        """Update both default and actual value of a config setting"""
        self.set_config(key, value)
        self._config[key]['default'] = self.get_config(key)

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
            1 CONFIG_LEVEL_ESSENTIAL: important setting which user has to consider
            2 CONFIG_LEVEL_BASIC: setting which many users may consider
            3 CONFIG_LEVEL_ADVANCED: setting relevant only to advanced users

        """

        import inspect
        import os

        caller = inspect.stack()[1]
        caller = os.path.splitext(os.path.basename(caller.filename))[0]
        logger.debug('Adding %i config items from %s' % (len(config), caller))
        remove = []
        for c, i in config.items():  # Check that provided config is conistent
            if c in self._config:
                if overwrite is False:
                    logger.debug(
                        '  Config item %s is already specified, not overwriting'
                        % c)
                    remove.append(c)
                else:
                    logger.debug('  Overwriting config item %s' % c)
            for p in ['type', 'description', 'level']:
                if p not in i:
                    raise ValueError(
                        '"%s" must be specified for config item %s' % (p, c))
            if i['level'] != CONFIG_LEVEL_ESSENTIAL and 'default' not in i:  #or i['default'] is None:
                raise ValueError(
                    'A default value must be provided for config item %s' % c)
            if i['type'] == 'enum':
                if 'enum' not in i or not isinstance(i['enum'], list):
                    raise ValueError(
                        '"enum" of type list must be provided for config item %s'
                        % (c))
            elif i['type'] in ['float', 'int']:
                for p in ['min', 'max', 'units']:
                    if p not in i:
                        raise ValueError(
                            '"%s" not provided for config item %s' % (p, c))
            elif i['type'] == 'bool':
                pass  # no check for bool
            else:
                raise ValueError(
                    'Config type "%s" (%s) is not defined. Valid options are: '
                    'float, int, enum, bool' % (i['type'], c))
            if 'default' in i:
                i['value'] = i['default']
        for r in remove:
            del config[r]
        self._config.update(config)

    def seed_elements(self,
                      lon,
                      lat,
                      time,
                      radius=0,
                      number=None,
                      radius_type='gaussian',
                      **kwargs):
        """Seed elements with given position(s), time and properties.

        Arguments:
            lon: scalar or array
                central longitude(s).
            lat: scalar or array
                central latitude(s).
            radius: scalar or array
                radius in meters around each lon-lat pair,
                within which particles will be randomly seeded.
            number: integer, total number of particles to be seeded
                If number is None, the number of elements is the
                length of lon/lat or time if these are arrays. Otherwise
                the number of elements are obtained from the config-default.
            time: datenum or list
                The time at which particles are seeded/released.
                If time is a list with two elements, elements are seeded
                continously from start/first to end/last time.
                If time is a list with more than two elements, the number of elements
                is equal to len(time) and are seeded as a time series.
            radius_type: string
                If 'gaussian' (default), the radius is the standard deviation in
                x-y-directions. If 'uniform', elements are spread evenly and
                always inside a circle with the given radius.
            kwargs:
                keyword arguments containing properties/attributes and
                values corresponding to the actual particle type (ElementType).
                These are forwarded to the ElementType class. All properties
                for which there are no default value must be specified.
        """

        if 'cone' in kwargs:
            raise ValueError(
                'Keyword *cone* for seed_elements is deprecated, use seed_cone() instead.'
            )

        if self.origin_marker is None:
            self.origin_marker = {}
        if 'origin_marker' in kwargs:
            origin_marker = kwargs['origin_marker']
        else:
            origin_marker = len(self.origin_marker)
        if 'origin_marker_name' in kwargs:
            origin_marker_name = kwargs['origin_marker_name']
            del kwargs['origin_marker_name']
        else:
            origin_marker_name = 'Seed %d' % len(self.origin_marker)
        if not 'origin_marker' in kwargs:
            kwargs['origin_marker'] = origin_marker
        if '_' in origin_marker_name:
            raise ValueError(
                'Underscore (_) not allowed in origin_marker_name')
        self.origin_marker[str(origin_marker)] = origin_marker_name.replace(
            ' ', '_')

        lon = np.atleast_1d(lon).ravel()
        lat = np.atleast_1d(lat).ravel()
        radius = np.atleast_1d(radius).ravel()
        time = np.atleast_1d(time)

        if lat.max() > 90 or lat.min() < -90:
            raise ValueError('Latitude must be between -90 and 90 degrees')

        if len(lon) != len(lat):
            raise ValueError('Lon and lat must have same lengths')

        if len(lon) > 1:
            if number is not None and number != len(lon):
                raise ValueError(
                    'Lon and lat have length %s, but number is %s' %
                    (len(lon), number))
            number = len(lon)
        else:
            if number is None:
                if len(time) > 2:
                    number = len(time)  # Interpreting as time series
                else:
                    number = self.get_config('seed:number')
            lon = lon * np.ones(number)
            lat = lat * np.ones(number)

        if len(time) != number and len(time) > 1:
            if len(time) == 2:  # start -> end
                td = (time[1] - time[0]) / (number - 1
                                            )  # timestep between points
                time = [time[0] + i * td for i in range(number)]
            else:
                raise ValueError(
                    'Time array has length %s, must be 1, 2 or %s' %
                    (len(time), number))

        # Add radius / perturbation
        if radius.max() > 0:
            geod = pyproj.Geod(ellps='WGS84')
            ones = np.ones(np.sum(number))
            if radius_type == 'gaussian':
                x = np.random.randn(np.sum(number)) * radius
                y = np.random.randn(np.sum(number)) * radius
                az = np.degrees(np.arctan2(x, y))
                dist = np.sqrt(x * x + y * y)
            elif radius_type == 'uniform':
                az = np.random.randn(np.sum(number)) * 360
                dist = np.sqrt(np.random.uniform(0, 1,
                                                 np.sum(number))) * radius
            lon, lat, az = geod.fwd(lon, lat, az, dist, radians=False)

        # If z is 'seafloor'
        if not 'z' in kwargs or kwargs['z'] is None:
            if 'seed:seafloor' in self._config:
                if self.get_config('seed:seafloor') is True:
                    kwargs['z'] = 'seafloor'
                    logger.debug('Seafloor is selected, neglecting z')
        if 'z' in kwargs and isinstance(kwargs['z'], str) \
                and kwargs['z'][0:8] == 'seafloor':
            # We need to fetch seafloor depth from reader
            seafloor_constant = self.get_config(
                'environment:constant:sea_floor_depth_below_sea_level')
            seafloor_fallback = self.get_config(
                'environment:fallback:sea_floor_depth_below_sea_level')
            if seafloor_constant is not None:
                env = {
                    'sea_floor_depth_below_sea_level':
                    np.array(seafloor_constant)
                }
            elif ('sea_floor_depth_below_sea_level'
                  in self.priority_list) or len(self._lazy_readers()):
                if not hasattr(self, 'time'):
                    self.time = time[0]
                env, env_profiles, missing = \
                    self.get_environment(['sea_floor_depth_below_sea_level'],
                                         time=time[0], lon=lon, lat=lat,
                                         z=0*lon, profiles=None)
            elif seafloor_fallback is not None:
                env = {
                    'sea_floor_depth_below_sea_level':
                    np.array(seafloor_fallback)
                }
            else:
                raise ValueError('A reader providing the variable '
                                 'sea_floor_depth_below_sea_level must be '
                                 'added before seeding elements at seafloor.')
            # Add M meters if given as 'seafloor+M'
            if len(kwargs['z']) > 8 and kwargs['z'][8] == '+':
                meters_above_seafloor = float(kwargs['z'][9::])
                logger.info('Seeding elements %f meters above seafloor' %
                            meters_above_seafloor)
            else:
                meters_above_seafloor = 0
            kwargs['z'] = \
                -env['sea_floor_depth_below_sea_level'].astype('float32') + meters_above_seafloor

        # Creating and scheduling elements
        elements = self.ElementType(lon=lon, lat=lat, **kwargs)
        time_array = np.array(time)
        self.schedule_elements(elements, time)

    def schedule_elements(self, elements, time):
        """Schedule elements to be seeded during runtime.

        Also assigns a unique ID to each particle, monotonically increasing."""

        # prepare time
        if isinstance(time, np.ndarray):
            time = list(time)
        if not isinstance(time, list):
            time = [time]
        if len(time) == 1 and len(elements) > 1:
            time = time * len(elements)

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
                logger.debug('Setting simulation start time to %s' %
                             str(min_time))
        else:
            self.start_time = min_time
            logger.debug('Setting simulation start time to %s' % str(min_time))

    @abstractproperty
    def ElementType(self):
        """Any trajectory model implementation must define an ElementType."""

    def split_along_elements(self, n='auto'):
        """
        Split init object along elements into separate simulations. Useful for parallel processing.

        Args:

            n: number of segments, defaults to 'auto': the number of cpus.
        """
        pass

    def copy(self):
        return copy.deepcopy(self)
