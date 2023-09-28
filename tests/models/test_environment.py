from opendrift.models.basemodel.environment import Environment
from opendrift.models.oceandrift import OceanDrift
from opendrift.config import Configurable


def test_add_readers(test_data_roms):
    c = Configurable()
    required_variables = {
        'x_sea_water_velocity': {
            'fallback': 0
        },
        'y_sea_water_velocity': {
            'fallback': 10
        },
    }
    env = Environment(required_variables, None, 1., c._config)
    env.add_reader(test_data_roms)
    env.finalize()

    e, p, k = env.get_environment(
        ['x_sea_water_velocity', 'y_sea_water_velocity'],
        test_data_roms.start_time, [10], [50], [0], [])
    assert e[0][0] == 0
    assert e[0][1] == 10


def test_fallback_values():
    o = OceanDrift()
    o.set_config('general:use_auto_landmask', False)
    assert o.required_variables == o.env.required_variables
    assert o._config == o.env._config
    assert len(o.env.fallback_values) == 0

    o.env.finalize()
    assert len(o.env.fallback_values) > 1

    assert o.env.fallback_values['x_sea_water_velocity'] == 0

    o.set_config('environment:fallback:x_sea_water_velocity', 10.)
    o.env.finalize()
    assert o.env.fallback_values['x_sea_water_velocity'] == 10.
