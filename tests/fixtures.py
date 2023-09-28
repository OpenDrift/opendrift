import pytest
import os
from pathlib import Path

@pytest.fixture
def test_data():
    import opendrift
    return os.path.abspath(
            os.path.join(os.path.dirname(opendrift.__file__),
                        '..', 'tests', 'test_data')) + os.path.sep

@pytest.fixture
def test_data_roms(test_data):
    from opendrift.readers import reader_ROMS_native
    return reader_ROMS_native.Reader(Path(test_data) / Path('2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc'))


@pytest.fixture
def show_plot(pytestconfig):
    return pytestconfig.getoption('plot')

