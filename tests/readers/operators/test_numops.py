import numpy as np
import pytest
from . import *
from opendrift.readers import reader_constant

def test_const_add():
    r = reader_constant.Reader({'wind_speed':5, 'wind_from_direction': 45, 'land_binary_mask': 0})
    r2 = r + 5
    print(r2)

    assert r2.variables == r.variables

    e, _ = r.get_variables_interpolated(variables = ['wind_speed'], lon=np.array([5.]), lat=np.array([60.]))
    w = e['wind_speed']
    print(w)

    e2, _ = r2.get_variables_interpolated(variables = ['wind_speed'], lon=np.array([5.]), lat=np.array([60.]))
    w2 = e2['wind_speed']
    print(w2)

    assert w2 == (w + 5)

