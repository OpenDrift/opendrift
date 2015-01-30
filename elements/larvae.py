import numpy as np

from elements import LagrangianArray

#########################################
# Generic Larvae class
#########################################


class Larvae(LagrangianArray):

    variables = LagrangianArray.add_variables(
        {'length':
            {'dtype': np.float32,
             'unit': 'mm'}})

#########################################
# Subclassing with specific Larva types
#########################################


class CodLarvae(Larvae):

    variables = Larvae.add_variables(
        {'CodLarvaeProperty1':
            {'dtype': np.float32}})


class HalibutLarvae(Larvae):

    variables = Larvae.add_variables(
        {'HalibutLarvaeProperty1':
            {'dtype': np.float32}})
