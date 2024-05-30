import logging
from typing import Dict, Union

logger = logging.getLogger(__name__)

CONFIG_LEVEL_ESSENTIAL = 1
CONFIG_LEVEL_BASIC = 2
CONFIG_LEVEL_ADVANCED = 3


class Configurable:
    _config: Union[Dict, None] = None

    def __init__(self):
        self._config = {}

    def get_config(self, key, default='raise'):
        if not key in self._config:
            if default == 'raise':
                raise ValueError('No config setting named %s' % key)
            else:
                return default
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
                elif i['type'] == 'str':
                    rang = f'min length {i["min_length"]}, max length {i["max_length"]}'
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
        if isinstance(value, dict):  # Recursive call with items in dictionary
            for subkey,subvalue in value.items():
                self.set_config(f'{key}:{subkey}', subvalue)
                logger.info(f'set_config(\'{key}:{subkey}\', {subvalue})')
            return
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
        elif i['type'] == 'str':
            if len(value) < i["min_length"] or len(value) > i["max_length"]:
                raise ValueError(f'String {key} length must be between {i["min_length"]} and {i["max_length"]} characters')
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

        type (string): 'float', 'int', 'str', 'bool' or 'enum'

        min, max (float/int/None): (only when type is 'float' or 'int')
            The minimum and maximum allowed values for this setting.
            May also be None if there are no upper/lowe limits.

        min_length, max_length (int): minimum and maximum length of string

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
            elif i['type'] == 'str':
                for p in ['min_length', 'max_length']:
                    if p not in i:
                        raise ValueError(
                            '"%s" not provided for config item %s' % (p, c))
            elif i['type'] == 'bool':
                pass  # no check for bool
            else:
                raise ValueError(
                    'Config type "%s" (%s) is not defined. Valid options are: '
                    'float, int, str, enum, bool' % (i['type'], c))
            if 'default' in i:
                i['value'] = i['default']
        for r in remove:
            del config[r]
        self._config.update(config)

