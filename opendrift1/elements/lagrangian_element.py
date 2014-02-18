
class LagrangianElement(object):

    variables = ['lon', 'lat', 'depth']


    @classmethod
    def update_properties(cls):
        print '\nUpdating properties of ' + cls.__name__ + ': ' + str(cls.variables)
