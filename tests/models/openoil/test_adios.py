import pytest
from opendrift.models.openoil import adios

@pytest.fixture
def aasgard():
    oils = adios.oils(1, 'AASGARD A 2003')
    f = oils[0].make_full()
    assert f.name == 'AASGARD A 2003'
    return f


def test_get_all():
    oils = adios.oils(300)
    print(oils[0])
    assert len(oils) == 300

    # make sure they are all unique
    ids = [o.id for o in oils]
    sids = set(ids)
    assert len(ids) == len(sids)

def test_get_full():
    thin = adios.oils(1)[0]
    full = thin.make_full()
    print(full)

def test_culled_density(aasgard):
    print(aasgard)
