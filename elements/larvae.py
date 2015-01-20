import numpy as np

from elements import LagrangianArray

#########################################
# Generic Larvae class
#########################################


class Larvae(LagrangianArray):

    parameters = LagrangianArray.add_parameters(
        {'length':
            {'dtype': np.float32,
             'unit': 'mm'}})

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
