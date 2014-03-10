import numpy as np

from elements import LagrangianArray

#########################################
# Generic Larvae class
#########################################


class Larvae(LagrangianArray):

    parameters = LagrangianArray.add_parameters(
        {'length':
            {'dtype': np.float32}})

    def update_properties(self):
        self.length = self.length*1.01  # General larvae grow by 1%


#########################################
# Subclassing with specific Larva types
#########################################

class CodLarvae(Larvae):

    parameters = Larvae.add_parameters(
        {'CodLarvaeProperty1':
            {'dtype': np.float32}})


class HalibutLarvae(Larvae):

    parameters = Larvae.add_parameters(
        {'HalibutLarvaeProperty1':
            {'dtype': np.float32}})

    def update_properties(self):
        self.length = self.length*1.02  # Halibut larvae grow by 2%
