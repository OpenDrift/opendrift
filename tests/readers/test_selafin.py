import os
import pytest

from . import *
from opendrift.readers.reader_telemac_selafin import Reader

has_telemac = False
try:
    from data_manip.formats.selafin import Selafin
    has_telemac = True
except ImportError:
    has_telemac = False

need_telemac = pytest.mark.skipif(has_telemac == False, reason = 'telemac pything bindings are not installed')

@pytest.fixture
def sel_3d(test_data):
    return os.path.join(test_data, 'Telemac_3D', 'r3d_tide_open_drift.slf')


# Lambert North
proj = "+proj=lcc +lat_1=49.50000000000001 +lat_0=49.50000000000001 +lon_0=0 \
        +k_0=0.999877341 +x_0=600000 +y_0=200000 +a=6378249.2 +b=6356515 \
        +units=m +no_defs"


@need_telemac
def test_open(sel_3d):
    r = Reader(sel_3d, proj4=proj)
    print(r)
