import numpy as np
import pytest
from . import *
from opendrift.readers import reader_constant

def test_const_add():
    r = reader_constant.Reader({'wind_speed':5, 'wind_from_direction': 45, 'land_binary_mask': 0})
    r2 = r + 5

