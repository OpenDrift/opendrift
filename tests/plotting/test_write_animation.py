import os
import pytest

from functools import lru_cache
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

@pytest.fixture
@lru_cache(maxsize = None)
def simulation():
    o = OceanDrift(loglevel=30)
    rn = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
    o.add_reader(rn)
    o.seed_elements(lon=4.8, lat=60.0, number=10, radius=1000,
                    time=rn.start_time)
    o.run(steps=5)
    return o

def test_write_mp4(simulation, tmpdir, benchmark):
    benchmark(simulation.animation, filename = tmpdir / 'test.mp4')

def test_write_gif(simulation, tmpdir, benchmark):
    benchmark(simulation.animation, filename = tmpdir / 'test.gif')
