from collections import OrderedDict
import logging

import numpy as np


class LagrangianArray(object):
    """A generic array-like class for Lagrangian particle tracking.

    A LagrangianArray is a generic class keeping the values of given
    properties ('variables') of a collection of particles at a given time.
    Values are stored as named attributes (similar to recarray) which are
    ndarrays (1D, vectors) with one value for each particle, or as scalars
    for values shared among all particles.

    This is an Abstract Base Class, meaning that only subclasses can be used.
    Subclasses will add specific variables for specific purposes (particle
    types, e.g. oil, fish eggs...) to the core variables described below.

    Attributes:
        variables: An OrderedDict where keys are names of the
            variables/properties of the current object. The values
            of the OrtderedDict are dictionaries with names such as
            'dtype', 'unit', 'standard_name' (CF), 'default' etc.
        All variable names will be added dynamically as attributes of
            the object after initialisation. These attributes will be
            numpy ndarrays of same length, or scalars. The core variables
            are:
        ID: an integer identifying each particle.
        status: 0 for active particles and a positive integer when deactivated
        lon: longitude (np.float32)
        lat: latitude (np.float32)
        depth: depth of the particle in m below sea surface.
    """

    variables = OrderedDict([
        ('ID', {'dtype': np.int32}),
        ('status', {'dtype': np.int32,
                    'default': 0}),
        ('lon', {'dtype': np.float32,
                 'units': 'degrees_east',
                 'standard_name': 'longitude',
                 'axis': 'X'}),
        ('lat', {'dtype': np.float32,
                 'units': 'degrees_north',
                 'standard_name': 'latitude',
                 'axis': 'Y'}),
        ('depth', {'dtype': np.float32,
                   'units': 'm',
                   'standard_name': 'depth',
                   'default': 0})])

    def __init__(self, **kwargs):
        """Initialises a LagrangianArray with given properties.

        Args:
            Keyword arguments (kwargs) with names corresponding to the
            OrderedDict 'variables' of the class, and corresponding values.
            The values must be ndarrays of equal length, or scalars.
            All (or none) variables must be given, unless a default value
            is specified in the OrderedDict 'variables'
            An empty object may be created by giving no input.
        """

        # Collect default values in separate dict, for easier access
        default_values = {variable: self.variables[variable]['dtype'](
                          self.variables[variable]['default'])
                          for variable in self.variables
                          if 'default' in self.variables[variable]}

        if len(kwargs) == 0:
            # Initialise an empty object (all variables have length 0)
            for var in self.variables:
                kwargs[var] = []

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
        for default_variable in default_values.keys():  # set default values
            setattr(self, default_variable, default_values[default_variable])
        for input_variable in kwargs.keys():  # override with input values
            setattr(self, input_variable, self.variables[input_variable]
                    ['dtype'](kwargs[input_variable]))

        # Store dtypes for all parameters in a common dtype object (for io)
        self.dtype = np.dtype([(var[0],var[1]['dtype'])
                               for var in self.variables.items()])

    @classmethod
    def add_variables(cls, new_variables):
        """Method used by subclasses to add specific properties/variables."""
        variables = cls.variables.copy()
        variables.update(new_variables)
        return variables

    def extend(self, other):
        """Add particles from another object. Scalar values are overwritten"""
        for var in self.variables:
            present_data = getattr(self, var)
            new_data = getattr(other, var)
            if isinstance(new_data, np.ndarray) and \
               isinstance(present_data, np.ndarray):  # Concatenating arrays
                setattr(self, var, np.concatenate((present_data, new_data)))
            else:
                setattr(self, var, new_data)  # NB: overwriting if scalar

    def move_elements(self, other, indices):
        """Remove elements with given indices, and append to another object."""
        # Move elements with given indices (boolean array)
        # to another LagrangianArray
        for var in self.variables:
            # Copy elements to other
            if isinstance(getattr(self, var), np.ndarray):  # Array
                setattr(other, var,
                        np.concatenate((getattr(other, var),
                                        getattr(self, var)[indices])))
                # Remove elements from self
                setattr(self, var, getattr(self, var)[~indices])
            else:
                setattr(other, var, getattr(self, var))  # Scalar
        if sum(indices) > 0:
            logging.debug('%s particles moved' % (sum(indices)))

    def __len__(self):
        return len(self.lat)

    def __repr__(self):
        outStr = ''
        for variable in self.variables.keys():
            outStr += variable + ': ' + str(getattr(self, variable)) + '\n'
        return outStr
