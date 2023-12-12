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

    # Oils in NOAA database with incomplete data, e.g. condensates without cuts. Should be updated
    to_be_fixed = ['LILLEFRIGG KONDENSAT 1996', 'SMORBUKK KONDENSAT 2003',
                   'GENERIC JET FUEL', 'GENERIC KEROSENE', 'MARINE DIESEL OIL, ESSO']
    to_be_fixed = []
    # Oils where name of sample is different from oil of reference
    # https://github.com/NOAA-ORR-ERD/noaa-oil-data/issues/1
    wrong_year = ['EKOFISK BLEND 2011', 'NORNE 1997', 'MIDGARD 1991']

    o = OpenOil(location = 'NORWAY')

    assert set(o.oiltypes).union(to_be_fixed).union(wrong_year) >= set(oils)
