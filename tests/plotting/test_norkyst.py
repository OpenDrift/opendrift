import pytest

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift


@pytest.mark.slow
@pytest.mark.mpl_image_compare
def test_fast_norkyst():
    o = OceanDrift(loglevel=30)
    rn = reader_netCDF_CF_generic.Reader(
        o.test_data_folder() +
        '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
    o.add_reader(rn)
    o.seed_elements(lon=4.8,
                    lat=60.0,
                    number=10,
                    radius=1000,
                    time=rn.start_time)
    o.run(steps=2)

    return o.plot(fast=True)[1]


@pytest.mark.slow
@pytest.mark.mpl_image_compare
def test_norkyst():
    o = OceanDrift(loglevel=30)
    rn = reader_netCDF_CF_generic.Reader(
        o.test_data_folder() +
        '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
    o.add_reader(rn)
    o.seed_elements(lon=4.8,
                    lat=60.0,
                    number=10,
                    radius=1000,
                    time=rn.start_time)
    o.run(steps=2)

    return o.plot(fast=False)[1]
