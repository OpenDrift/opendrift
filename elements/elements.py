from abc import abstractmethod, ABCMeta
from collections import OrderedDict

import numpy as np


class LagrangianArray(object):
    __metaclass__ = ABCMeta

    variables = OrderedDict([
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
        default_values = {variable: self.variables[variable]['dtype'](
                          self.variables[variable]['default'])
                          for variable in self.variables
                          if 'default' in self.variables[variable]}

        # Check for missing arguments
        missing_args = set(self.variables.keys()) - \
            set(kwargs.keys()) - set(default_values.keys())
        if missing_args:
            raise TypeError('Missing arguments: %s' % str(list(missing_args)))

        # Check for redundant arguments
        redundant_args = set(kwargs.keys() + default_values.keys()) - set(
            self.variables.keys())
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
        # NB: should check that variable names are not reserved!
        for default_variable in default_values.keys():  # set default values
            setattr(self, default_variable, default_values[default_variable])
        for input_variable in kwargs.keys():  # override with input values
            # presently as scalars, but should be arrays?
            setattr(self, input_variable, self.variables[input_variable]
                    ['dtype'](kwargs[input_variable]))

    @classmethod
    def add_variables(cls, new_variables):
        variables = cls.variables.copy()
        variables.update(new_variables)
        return variables

    def __len__(self):
        return len(self.lat)

    def show_data(self):
        for variable in self.variables.keys():
            print variable + ': ' + str(getattr(self, variable))
