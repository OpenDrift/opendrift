import pytest
from opendrift.models.openoil import adios


@pytest.fixture
def aasgard():
    oils = adios.oils(1, 'AASGARD A 2003')
    f = oils[0].make_full()
    assert f.name == 'AASGARD A 2003'
    return f


def test_culled_density(aasgard):
    print(aasgard)
