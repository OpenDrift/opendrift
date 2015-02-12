import numpy as np

from elements import LagrangianArray

#########################################
# Generic Larvae class
#########################################


class Larvae(LagrangianArray):
    """Extending LagrangianArray with variables relevant for (marine) larvae."""

    variables = LagrangianArray.add_variables(
        {'length':
            {'dtype': np.float32,
             'unit': 'mm'}})

#########################################
# Subclassing with specific Larva types
#########################################


class CodLarvae(Larvae):
    """Extending Larvae with variables relevant for cod larvae."""

    variables = Larvae.add_variables(
        {'CodLarvaeProperty1':
            {'dtype': np.float32}})


class HalibutLarvae(Larvae):
    """Extending Larvae with variables relevant for halibut larvae."""

    variables = Larvae.add_variables(
        {'HalibutLarvaeProperty1':
            {'dtype': np.float32}})
