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
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
# 
# Copyright 2015, Knut-Frode Dagestad, MET Norway

from opendrift.elements import LagrangianArray


class PassiveTracer(LagrangianArray):
    """Basic implementation of LagrangianArray with no additional properties.

    Contains only the properties of the abstract class LagrangianArray,
    i.e. position (lon, lat, z) and ID.
    May be used for passive tracer calculations when no properties are needed.
    """

    pass
