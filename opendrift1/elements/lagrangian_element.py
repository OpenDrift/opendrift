
class LagrangianElement(object):

    # Variables are updated for all elements at each time step
    variables = [('lon', 'float32'),
                 ('lat', 'float32'),
                 ('depth', 'float32')]

    # Constants are constant for each element
    constants = [('particleID', 'int32')] # particle ID, integer


    @classmethod
    def update_properties(cls):
        print '\nUpdating properties of ' + cls.__name__ + \
                ': ' + str(cls.variables)
