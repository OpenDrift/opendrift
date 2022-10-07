import pytest

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.basemodel import *
from opendrift.models.basemodel.environment import Environment
from opendrift.models.basemodel.config import Configurable

def test_add_reader(test_data):
    reader_arome = reader_netCDF_CF_generic.Reader(test_data / '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
    env = Environment([])
    env.add_reader(reader_arome)

@pytest.mark.slow
def test_auto_landmask():
    env = Environment({ 'land_binary_mask': {'fallback': None}, })
    c = Configurable()
    c._add_config({
        # type, default, min, max, enum, important, value, units, description
        'general:use_auto_landmask': {
            'type':
            'bool',
            'default':
            True,
            'description':
            'A built-in GSHHG global landmask is used if True, '
            'otherwise landmask is taken from reader or fallback value.',
            'level':
            CONFIG_LEVEL_ADVANCED
        }})
    env2 = env.finalize(c)
    print(env2.list_environment_variables())

    assert 'land_binary_mask' in env2.list_environment_variables()

def test_missing_variable():
    env = Environment({ 'x_wind': {'fallback': None}, })
    c = Configurable()
    c._add_config({
        # type, default, min, max, enum, important, value, units, description
        'general:use_auto_landmask': {
            'type':
            'bool',
            'default':
            True,
            'description':
            'A built-in GSHHG global landmask is used if True, '
            'otherwise landmask is taken from reader or fallback value.',
            'level':
            CONFIG_LEVEL_ADVANCED
        }})
    with pytest.raises(ValueError):
        env = env.finalize(c)

