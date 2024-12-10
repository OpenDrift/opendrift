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
    # previous oils
    previous_oils_file = Path(test_data) / 'generated' / 'oil_list_norway.txt'
    with open(previous_oils_file, 'r') as fd:
        previous_oils = [o.strip() for o in fd.readlines()]

    o = OpenOil(location = 'NORWAY')

    new_oils = set(o.oiltypes) - set(previous_oils)
    removed_oils = set(previous_oils) - set(o.oiltypes)

    if len(new_oils) > 0:
        print('\nWarning, new oils:')
        for ot in new_oils:
            print(ot)

    if len(removed_oils) > 0:
        print('\nWarning, removed oils:')
        for ot in removed_oils:
            print(ot)

    if False:  # Change to True when updating list of oils
        print(f'Writing present list of oiltypes to {previous_oils_file}')
        with open(previous_oils_file, 'w') as f:
            for oil in o.oiltypes:
                f.write(f"{oil}\n")
    
    assert len(new_oils) == 0
    assert len(removed_oils) == 0
