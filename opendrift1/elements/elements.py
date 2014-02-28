import numpy as np

class LagrangianArray(object):

    # Variables are updated for all elements at each time step
    # name, data type[, default value]
    variables = [('lon', 'float32'),
                 ('lat', 'float32'),
                 ('depth', 'float32', 0)]

    # Constants are constant for each element
    constants = [('particleID', 'int32')] # particle ID, integer


    def __init__(self, **kwargs):

        # Reformat variables information
        # Allowing simplest possible specification (list) in modules
        self.variable_names = [var[0] for var in self.variables] 
        self.variable_dtypes = [var[1] for var in self.variables] 
        self.variable_defaults = {
            var[0]: var[2] for var in self.variables if len(var) == 3}
        # Check for missing arguments
        missing_args = set(self.variable_names) - \
            set(kwargs.keys()) - set(self.variable_defaults.keys())
        if missing_args:
            raise TypeError('Missing argument[s]: %s' % str(list(missing_args)))

        # Check for redundant arguments
        redundant_args = set(kwargs.keys() + self.variable_defaults.keys()) - set(self.variable_names)
        if redundant_args:
            raise TypeError('Redundant argument[s]: %s' % str(list(redundant_args)))

        # Check that input arrays have same length
        array_lengths = set([len(array) for array in kwargs.values()])
        if len(array_lengths - {1}) > 1:
            raise TypeError(
                'Input arrays must have same length. Lengths given: ' \
                + str(list(array_lengths)))
        

    def update_properties(self):
        print '\nUpdating properties:' + \
                ': ' + str(self.variables)
