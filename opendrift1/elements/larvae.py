import numpy as np

from elements import LagrangianArray

#########################################
# Generic Larvae class
#########################################

class Larvae(LagrangianArray):

    variables = LagrangianArray.variables + [('length', np.float32)]

    def update_properties(cls):
        super(Larvae, cls).update_properties()
        print '...using specialised function for Larvae class'


#########################################
# Subclassing with specific Larva types 
#########################################

class CodLarvae(Larvae):

    variables = Larvae.variables + [
                                ('CodLarvaeProperty1', np.float32)]

    def update_properties(self):
        super(Larvae, self).update_properties()
        print '...using specialised function for CodLarvae class'



class HalibutLarvae(Larvae):

    variables = Larvae.variables + [
                    ('HalibutLarvaeProperty1', np.float32)]

    def update_properties(self):
        super(Larvae, self).update_properties()
        print '...using specialised function for HalibutLarvae class'


