# This file is part of OpenDrift.
# 
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
# 
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
# 
# Copyright 2015, Knut-Frode Dagestad, MET Norway

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
