import numpy as np
import pytest
from . import *
from opendrift.readers import reader_global_landmask

def test_global_setup(benchmark):
    benchmark(reader_global_landmask.Reader)

@pytest.mark.slow
def test_performance_global(benchmark):
    print("setting up global landmask")
    reader_global = reader_global_landmask.Reader(
            extent = [18.64, 69.537, 19.37, 69.81])


    x = np.linspace(18.641, 19.369, 100)
    y = np.linspace(69.538, 69.80, 100)

    xx, yy = np.meshgrid(x,y)
    xx = xx.ravel()
    yy = yy.ravel()

    print ("points:", len(xx))

    # warmup
    reader_global.__on_land__(xx, yy)

    print ("masking against cartopy")
    benchmark(reader_global.__on_land__, xx,yy)

