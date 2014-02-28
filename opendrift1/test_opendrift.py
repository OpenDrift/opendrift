#!/usr/bin/env python

import numpy as np

from elements import LagrangianArray
from elements.oil import Oil
from elements.larvae import Larvae, CodLarvae, HalibutLarvae


o = Oil(lon=[32, 3, 8], lat=np.array([22, 3, 4]), depth=[44], massOil=[99])
#o = Oil(lat=[22, 3, 4], depth=[44], massEvaporated=[99])
#o.update_properties()
#print o
#
c = CodLarvae(lon=[32, 3, 8], lat=[22, 3, 4], depth=[44], CodLarvaeProperty1=[5], length=[10])
#c.update_properties()
#print c
#print o.variables['depth']
