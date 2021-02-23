import pytest
from . import *

@pytest.mark.skip("missing suitable data")
def test_load_grb_file(test_data):
    from opendrift.readers.reader_grib import Reader
    waver = Reader(test_data + 'norkyst800m_weatherapi_nordland.grb')
    print("reader:", waver)

