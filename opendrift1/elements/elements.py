import numpy as np

class LagrangianArray(object):

    parameters = {
        'lon': {'dtype': np.float32, 
                'unit': 'degrees',
                'short_name': 'longitude'},
        'lat': {'dtype': np.float32,
                'unit': 'degrees',
                'short_name': 'latitude'},
        'depth': {'dtype': np.float32,
                  'unit': 'm',
                  'short_name': 'depth',
                  'default': 0},
        'particleID': {'dtype': np.int32,
                       'time_independent': True}}



    def __init__(self, **kwargs):
        # Reformat variables information
        # Allowing simplest possible specification (list) in modules
        #self.variable_dtypes = [var[1] for var in self.variables] 
        #self.variable_defaults = {
        #    var[0]: var[1](var[2]) for var in self.variables if len(var) == 3}
        print self.parameters
        print self.parameters.keys()
        print type(self.parameters)
        print type(self.parameters.keys())
        #sys.exit('stop')
        default_values = {par.key: par['dtype'](par['default'])
                            for par in self.parameters
                            if par.has_key('default')}

        print default_values
        return
        

        # Check for missing arguments
        missing_args = set(self.parameters.keys()) - \
        set(kwargs.keys()) - set(self.variable_defaults.keys())
        if missing_args:
            raise TypeError('Missing argument[s]: %s' % str(list(missing_args)))

        # Check for redundant arguments
        redundant_args = set(kwargs.keys() + self.variable_defaults.keys()) - set(self.variable_names)
        if redundant_args:
            raise TypeError('Redundant argument[s]: %s' % str(list(redundant_args)))

        # Check that input arrays have same length
        array_lengths =  [1]*len(kwargs)
        for i, input_value in enumerate(kwargs.values()):
            try:
                array_lengths[i] = len(array)
            except:
                array_lengths[i] = 1 # scalar is given
        if len(set(array_lengths) - {1}) > 1:
            raise TypeError(
                'Input arrays must have same length. Lengths given: ' \
                + str(array_lengths))

        # Store input arrays
        # NB: should check that dynamic names are not reserved!
        for default_name, default_value in self.variable_defaults.iteritems(): # set default values
            setattr(self, default_name, default_value)
        for i, key in enumerate(kwargs.keys()): # override with input values
            # presently as scalars, but should be arrays?
            setattr(self, key, self.variable_dtypes[i](kwargs[key]))
        
    def __len__(self):
        return len(self.lat)

    def show_data(self):
        for varname in self.variable_names:
            print varname + ': ' + str(getattr(self, varname))


    def update_properties(self):
        print '\nUpdating properties:' + \
                ': ' + str(self.variables)
