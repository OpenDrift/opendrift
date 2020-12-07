import numpy as np
import pytest
import matplotlib.pyplot as plt
from tests import *
from opendrift.readers.reader_grib2 import Reader as Grib
from opendrift.models.oceandrift import OceanDrift

caps_u10 = 'grb/u_vel.grb2'
caps_v10 = 'grb/v_vel.grb2'

def test_open():
    r = Grib(caps_u10)
    print(r)
    # r.plot_mesh()
    # plt.show()

