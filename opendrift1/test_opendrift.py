#!/usr/bin/env python

from elements.lagrangian_element import LagrangianElement
from elements.oil import OilElement
from elements.larvae import LarvaeElement
from state import State

LagrangianElement.update_properties()

OilElement.update_properties()

#LarvaeElement.update_properties()

s = State(OilElement, x=32, y=22)
