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

from collections import OrderedDict
import numpy as np

class LagrangianArray:
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
            of the OrderedDict are dictionaries with names such as
            'dtype', 'unit', 'standard_name' (CF), 'default' etc.
            All variable names will be added dynamically as attributes of
            the object after initialisation. These attributes will be
            numpy ndarrays of same length, or scalars. The core variables
            are:

                - ID: an integer identifying each particle.
                - status: 0 for active particles and a positive integer when deactivated
                - lon: longitude (np.float32)
                - lat: latitude (np.float32)
                - z: vertical position of the particle in m, positive upwards (above sea surface)
    """

    variables = OrderedDict([
        ('ID', {'dtype': np.int32,  # Unique numerical identifier
                'seed': False,
                'default': -1}),  # ID to be assigned by application
        ('status', {'dtype': np.int32,  # Status categories
                    'seed': False,
                    'default': 0}),
        ('moving', {'dtype': np.int32,  # Set to 0 for elements which are frosen
                    'seed': False,
                    'default': 1}),
        ('age_seconds', {'dtype': np.float32,
                         'units': 's',
                         'seed': False,
                         'default': 0}),
        ('origin_marker', {'dtype': np.int32,
                           'unit': '',
            'description': 'An integer kept constant during the simulation. Different values may be used for different seedings, to separate elements during analysis. With GUI, only a single seeding is possible.',
                           'default': 0}),
        ('lon', {'dtype': np.float32,
                 'units': 'degrees_east',
                 'standard_name': 'longitude',
                 'long_name': 'longitude',
                 'seed': False,
                 'axis': 'X'}),
        ('lat', {'dtype': np.float32,
                 'units': 'degrees_north',
                 'standard_name': 'latitude',
                 'long_name': 'latitude',
                 'seed': False,
                 'axis': 'Y'}),
        ('z', {'dtype': np.float32,
                   'units': 'm',
                   'standard_name': 'z',
                   'long_name': 'vertical position',
                   'axis': 'Z',
                   'positive': 'up',
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
        redundant_args = set(list(kwargs.keys()) +
            list(default_values.keys())) - set((self.variables.keys()))
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
        self.dtype = np.dtype([(var[0], var[1]['dtype'])
                               for var in self.variables.items()])

        # Status must always be array
        if not type(self.status) == np.ndarray:
            self.status = self.status*np.ones(self.lon.shape)

    @classmethod
    def add_variables(cls, new_variables):
        """Method used by subclasses to add specific properties/variables."""
        variables = cls.variables.copy()
        variables.update(new_variables)
        return variables

    def extend(self, other):
        """Add elements from another object."""
        len_self = len(self)
        len_other = len(other)
        for var in self.variables:
            present_data = getattr(self, var)
            new_data = getattr(other, var)

            # If both arrays have an identical scalar, it remains a scalar
            if (not isinstance(new_data, np.ndarray) and
                not isinstance(present_data, np.ndarray) and
                present_data == new_data):
                continue

            else:  # Otherwise we create arrays and concatenate
                if not hasattr(present_data, '__len__'):
                    present_data = present_data*np.ones(len_self)
                if not hasattr(new_data, '__len__'):
                    new_data = new_data*np.ones(len_other)
                setattr(self, var, np.concatenate((present_data,
                                                   new_data)))

    def move_elements(self, other, indices):
        """Remove elements with given indices, and append to another object.
        NB: indices is boolean array, not real indices!"""

        # Move elements with given indices (boolean array)
        # to another LagrangianArray
        # NB: scalars and 1D arrays are converted to ndarrays and concatenated
        self_len = len(self)
        other_len = len(other)
        for var in self.variables:
            self_var = getattr(self, var)
            other_var = getattr(other, var)
            if (not isinstance(self_var, np.ndarray) and
                not isinstance(other_var, np.ndarray)) and \
                    (other_var == self_var):
                    if np.sum(indices) == len(self):
                        setattr(self, var, [])  # Empty if all elements moved
                    continue  # Equal scalars - we do nothing

            # Copy elements to other
            self_var = np.atleast_1d(self_var)
            other_var = np.atleast_1d(other_var)
            if len(self_var) < self_len:  # Convert scalar to array
                self_var = self_var*np.ones(self_len)
            if len(other_var) < other_len:  # Convert scalar to aray
                other_var = other_var*np.ones(other_len)
            if len(self_var) > 0:
                setattr(other, var, np.concatenate((other_var,
                                                    self_var[indices])))
            else:
                setattr(other, var, self_var[indices])
            setattr(self, var, self_var[~indices])  # Remove from self

            #if isinstance(self_var, np.ndarray) or\
            #    isinstance(other_var, np.ndarray):  # Array
            #    setattr(other, var,
            #            np.concatenate((getattr(other, var),
            #                            getattr(self, var)[indices])))
            #    # Remove elements from self
            #    setattr(self, var, getattr(self, var)[~indices])
            #elif isinstance(getattr(other, var), np.ndarray):
            #    setattr(other, var,
            #            np.concatenate((getattr(other, var),
            #                            np.atleast_1d(getattr(self, var)))))
            #else:
            #    setattr(other, var, getattr(self, var))  # Scalar

    def __len__(self):
        length = 0
        for var in self.variables:
            length = np.maximum(length, len(np.atleast_1d(getattr(self, var))))
        return length

    def __repr__(self):
        outStr = ''
        for variable in self.variables.keys():
            outStr += variable + ': ' + str(getattr(self, variable)) + '\n'
        return outStr
