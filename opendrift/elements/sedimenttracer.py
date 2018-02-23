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
from  buoyanttracer import BuoyantTracer
import numpy as np

class SedimentTracer(BuoyantTracer):
    """Extending BuoyantTracer (which is same as PassiveTracer+settling velocity) to a sediment element
       which has additional properties such as critical_shear_stress for resuspension 
    """

    variables = BuoyantTracer.add_variables([
        ('critical_shear_stress', {'dtype': np.float32,
                               'units': 'm/s',
                               'default': 0.}),

        # ('Other??', {'dtype': np.float32,
        #                        'unit': '%',
        #                        'default': 0.0}),

        ])