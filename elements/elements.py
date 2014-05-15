from abc import abstractmethod, ABCMeta
from collections import OrderedDict

import numpy as np


class LagrangianArray(object):
    __metaclass__ = ABCMeta

    parameters = OrderedDict([
        ('lon', {'dtype': np.float32,
                 'unit': 'degrees',
                 'standard_name': 'longitude'}),
        ('lat', {'dtype': np.float32,
                 'unit': 'degrees',
                 'standard_name': 'latitude'}),
        ('depth', {'dtype': np.float32,
                   'unit': 'm',
                   'standard_name': 'depth',
                   'default': 0})])

    def __init__(self, **kwargs):

        # Collect default values in separate dict
        default_values = {parameter: self.parameters[parameter]['dtype'](
                          self.parameters[parameter]['default'])
                          for parameter in self.parameters
                          if 'default' in self.parameters[parameter]}

        # Check for missing arguments
        missing_args = set(self.parameters.keys()) - \
            set(kwargs.keys()) - set(default_values.keys())
        if missing_args:
            raise TypeError('Missing arguments: %s' % str(list(missing_args)))

        # Check for redundant arguments
        redundant_args = set(kwargs.keys() + default_values.keys()) - set(
            self.parameters.keys())
        if redundant_args:
            raise TypeError('Redundant arguments: %s' %
                            str(list(redundant_args)))

        # Check that input arrays have same length
        array_lengths = [1]*len(kwargs)
        for i, input_value in enumerate(kwargs.values()):
            try:
                array_lengths[i] = len(input_value)
            except:
                array_lengths[i] = 1  # scalar is given
        if len(set(array_lengths) - {1}) > 1:
            raise TypeError(
                'Input arrays must have same length. Lengths given: '
                + str(array_lengths))

        # Store input arrays
        # NB: should check that parameter names are not reserved!
        for default_parameter in default_values.keys():  # set default values
            setattr(self, default_parameter, default_values[default_parameter])
        for input_parameter in kwargs.keys():  # override with input values
            # presently as scalars, but should be arrays?
            setattr(self, input_parameter, self.parameters[input_parameter]
                    ['dtype'](kwargs[input_parameter]))

    @classmethod
    def add_parameters(cls, new_parameters):
        parameters = cls.parameters.copy()
        parameters.update(new_parameters)
        return parameters

    def __len__(self):
        return len(self.lat)

    def show_data(self):
        for parameter in self.parameters.keys():
            print parameter + ': ' + str(getattr(self, parameter))

    @abstractmethod
    def update_properties(self, environment):
        print '\nUpdating properties: ' + \
            str(self.parameters.keys())

    @abstractmethod
    def update_position(self, readers):
        print '\nUpdating position.'
