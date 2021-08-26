import pytest
from opendrift.models.openoil import adios


def test_get_300():
    oils = adios.oils(300)
    print(oils[0])
    assert len(oils) == 300

    # make sure they are all unique
    ids = [o.id for o in oils]
    sids = set(ids)
    assert len(ids) == len(sids)

@pytest.mark.skip
def test_get_all():
    oils = adios.oils(None)
    print(oils[0])
    assert len(oils) > 1800

    # make sure they are all unique
    ids = [o.id for o in oils]
    sids = set(ids)
    assert len(ids) == len(sids)

def test_get_full():
    thin = adios.oils(1)[0]
    full = thin.make_full()
    print(full)

def test_bunker_c_1987():
    o = adios.oils(query='Bunker C [1987]')
    print([oo.name for oo in o])
    assert len(o) > 1
    o = adios.find_full_oil_from_name('Bunker C [1987]')
    assert o.valid()

