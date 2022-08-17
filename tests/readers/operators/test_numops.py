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

    e2, _ = r2.get_variables_interpolated(variables = ['wind_speed'], lon=np.array([5.]), lat=np.array([60.]))
    w2 = e2['wind_speed']

    assert w2 == (w + 5)

def test_reader_add():
    r1 = reader_constant.Reader({'wind_speed':5, 'wind_from_direction': 45, 'land_binary_mask': 0})
    r2 = reader_constant.Reader({'wind_speed':7, 'wind_from_direction': 30, 'land_binary_mask': 0})

    r = r1.filter_vars(['wind_speed']) + r2
    e, _ = r.get_variables_interpolated(variables = ['wind_speed'], lon=np.array([5.]), lat=np.array([60.]))
    w = e['wind_speed']
    assert w == [(5+7)]

    r = r1.exclude_vars(['wind_from_direction']) + r2

    with pytest.raises(AssertionError):
        e, _ = r.get_variables_interpolated(variables = ['wind_from_direction'], lon=np.array([5.]), lat=np.array([60.]))

def test_reader_div():
    r1 = reader_constant.Reader({'wind_speed':5, 'wind_from_direction': 45, 'land_binary_mask': 0})
    r2 = reader_constant.Reader({'wind_speed':7, 'wind_from_direction': 30, 'land_binary_mask': 0})

    r = (r1.filter_vars(['wind_from_direction']) + r2) / 2
    e, _ = r.get_variables_interpolated(variables = ['wind_from_direction'], lon=np.array([5.]), lat=np.array([60.]))
    w = e['wind_from_direction']
    assert w == [(75/2)]


def filter_vars():
    r1 = reader_constant.Reader({'wind_speed':5, 'wind_from_direction': 45, 'land_binary_mask': 0})
    assert r1.variables == [ 'wind_speed', 'wind_from_direction', 'land_binary_mask' ]

    r2 = r1.exclude_vars(['wind_from_direction'])
    assert r2.variables == [ 'wind_speed', 'land_binary_mask' ]

