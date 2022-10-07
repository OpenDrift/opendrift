import pytest
import os
from pathlib import Path

@pytest.fixture
def test_data():
    import opendrift
    return Path(os.path.abspath(
            os.path.join(os.path.dirname(opendrift.__file__),
                        '..', 'tests', 'test_data')))

@pytest.fixture
def show_plot(pytestconfig):
    return pytestconfig.getoption('plot')

