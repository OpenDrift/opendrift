import pytest
from pathlib import Path
from opendrift.models.openoil import OpenOil

@pytest.mark.xfail
def test_get_openoil_list(test_data):
    oils = Path(test_data) / 'generated' / 'oil_list_full.txt'
    with open(oils, 'r') as fd:
        oils = [o.strip() for o in fd.readlines()]

    o = OpenOil()

    assert set(o.oiltypes) >= set(oils)

def test_get_openoil_list_norway(test_data):
    oils = Path(test_data) / 'generated' / 'oil_list_norway.txt'
    with open(oils, 'r') as fd:
        oils = [o.strip() for o in fd.readlines()]

    o = OpenOil(location = 'NORWAY')

    print(set(oils) - set(o.oiltypes))

    assert set(o.oiltypes) >= set(oils)
