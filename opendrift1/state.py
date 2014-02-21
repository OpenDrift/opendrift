import numpy as np

class State():

    def __init__(self, LagrangianElementClass, **kwargs):

        self.dtypes = LagrangianElementClass.variables

        for key, value in kwargs.iteritems():
            print key, value


        self.variables = np.ndarray()
        pass

    def __repr__(self):
        return str(self.dtypes)
