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

# from elements import LagrangianArray
from  passivetracer import PassiveTracer
import numpy as np

class BuoyantTracer(PassiveTracer):
    """Extending PassiveTracer (which is same as LagrangianArray) for elements moving in 3 dimensions
       The BuoyantTracer may be buoyant and/or subject to vertical mixing, and has some additional properties (as used in oceandrift3D.py)
       Buoyant behaviour is described by terminal velocity
       
       terminal_velocity>0 particle moves up towards the surface
       terminal_velocity<0 particle moves down towards the seabed
    """

    variables = PassiveTracer.add_variables([
        ('terminal_velocity', {'dtype': np.float32,
                               'units': 'm/s',
                               'default': 0.}),

        ('wind_drift_factor', {'dtype': np.float32,
                               'unit': '%',
                               'default': 0.0}),

        ('age_seconds', {'dtype': np.float32,
                         'units': 's',
                         'default': 0})

        ])