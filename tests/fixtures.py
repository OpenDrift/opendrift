import pytest
import os

@pytest.fixture
def test_data():
    import opendrift
    return os.path.abspath(
            os.path.join(os.path.dirname(opendrift.__file__),
                        '..', 'tests', 'test_data')) + os.path.sep

@pytest.fixture
def show_plot(pytestconfig):
    return pytestconfig.getoption('plot')

