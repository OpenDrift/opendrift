import numpy as np

from elements import LagrangianArray

# Modify Class name and variable name/properties
# below to make a new element class


class PassiveTracer(LagrangianArray):

    def update_properties(self):
        # overload from superclass
        super(PassiveTracer, self).update_properties()
        pass
