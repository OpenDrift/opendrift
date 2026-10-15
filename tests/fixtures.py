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
def test_data_fake_ig(test_data):
    """
    Return reader_netCDF_CF_generic.Reader containing fake wind data
    with same coverage as the DMI HARMONIE/IG model (Greenland)
    """
    from opendrift.readers import reader_netCDF_CF_generic
    return reader_netCDF_CF_generic.Reader(Path(test_data) / 'harmonie_ig_coverage' / 'fake_dmi_ig_wind.nc')

@pytest.fixture
def test_data_fake_wind_ob_tran(test_data):
    """
    Return reader_netCDF_CF_generic.Reader containing fake wind data
    with same coverage as the old DMI NEA model (using rotated pole / ob_tran projection)
    """
    from opendrift.readers import reader_netCDF_CF_generic
    return reader_netCDF_CF_generic.Reader(Path(test_data) / 'ob_tran_coverage' / 'fake_wind_ob_tran.nc')

@pytest.fixture
def show_plot(pytestconfig):
    return pytestconfig.getoption('plot')

