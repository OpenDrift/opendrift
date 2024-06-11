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
    env = Environment(required_variables, c._config)
    env.add_reader(test_data_roms)
    env.finalize()

    e, p, k = env.get_environment(
        ['x_sea_water_velocity', 'y_sea_water_velocity'],
        test_data_roms.start_time, [10], [50], [0], [])
    assert e[0][0] == 0
    assert e[0][1] == 10
