from opendrift.models.basemodel.environment import Environment
from opendrift.config import Configurable

def test_add_readers(test_data_roms):
    c = Configurable()
    required_variables = {
        'x_sea_water_velocity': {'fallback': 0}, }
    env = Environment(required_variables, None, 1., c._config)
    env.add_reader(test_data_roms)
    env.finalize()

    e, p = env.get_environment(['x_sea_water_velocity'], test_data_roms.start_time, [10], [50], [0], [])

