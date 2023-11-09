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

    print(oils, 'OILS')

    o = OpenOil(location = 'NORWAY')
    print(o.oiltypes, 'OILTYPES')
    print(len(o.oiltypes))
    print(len(oils))
    print(set(o.oiltypes)-set(oils))

    assert set(o.oiltypes) >= set(oils)
