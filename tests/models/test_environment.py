from datetime import datetime, timedelta
import pytest
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

def test_previous():
    o = OceanDrift()  # With coastline_action = none
    o.set_config('general:coastline_action', 'none')
    o.set_config('environment:constant:land_binary_mask', 0)
    o.set_config('environment:constant:x_sea_water_velocity', 1)
    o.seed_elements(lon=3, lat=60, time=datetime.now())
    o.run(steps=1)
    assert o.elements.lon == pytest.approx(3.0645, .001)
    assert o.elements_previous is None
    assert o.environment_previous is None

    o = OceanDrift()  # With coastline_action = previous
    o.set_config('general:coastline_action', 'previous')
    o.set_config('environment:constant:land_binary_mask', 0)
    o.set_config('environment:constant:x_sea_water_velocity', 1)
    o.seed_elements(lon=3, lat=60, time=datetime.now())
    o.run(steps=1)
    assert o.elements.lon == pytest.approx(3.0645, .001)
    assert o.elements_previous.lon == pytest.approx(3.0, .001)
    assert o.elements_previous is not None
    assert o.environment_previous is None
 
def test_skip_env_variable():
    o = OceanDrift()
    o.set_config('drift:vertical_mixing', True)  # Diffusivity shall be included
    o.set_config('environment:constant:land_binary_mask', 0)
    o.seed_elements(lon=3, lat=60, time=datetime.now())
    o.run(steps=1)
    assert 'ocean_vertical_diffusivity' in o.required_variables

    o = OceanDrift()
    o.set_config('drift:vertical_mixing', False)  # Diffusivity shall be skipped
    o.set_config('environment:constant:land_binary_mask', 0)
    o.seed_elements(lon=3, lat=60, time=datetime.now())
    o.run(steps=1)
    assert 'ocean_vertical_diffusivity' not in o.required_variables
