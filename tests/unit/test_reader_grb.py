import pytest
from . import *

def test_load_grb_file(test_data):
    from opendrift.readers.reader_grib import Reader
    waver = Reader(test_data + 'mywavewam4_weatherapi_n-northsea.grb')
    print("reader:", waver)

