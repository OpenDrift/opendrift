import pytest

from datetime import datetime
from opendrift import test_data_folder as tdf
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift

snooze_time = datetime(2026, 7, 1)
snooze_graphical = datetime.now() < snooze_time
snooze_message = f'Snoozing graphical tests until {snooze_time}'

@pytest.mark.skipif(snooze_graphical is True, reason=snooze_message)
@pytest.mark.slow
@pytest.mark.mpl_image_compare
def test_fast_norkyst():
    o = OceanDrift(loglevel=30)
    rn = reader_netCDF_CF_generic.Reader(
        tdf +
        '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
    o.add_reader(rn)
    o.seed_elements(lon=4.8,
                    lat=60.0,
                    number=10,
                    radius=1000,
                    time=rn.start_time)
    o.run(steps=2)

    return o.plot(fast=True, buffer=.2, land_zorder=0)[1]


@pytest.mark.skipif(snooze_graphical is True, reason=snooze_message)
@pytest.mark.slow
@pytest.mark.mpl_image_compare
def test_norkyst():
    o = OceanDrift(loglevel=30)
    rn = reader_netCDF_CF_generic.Reader(
        tdf +
        '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
    o.add_reader(rn)
    o.seed_elements(lon=4.8,
                    lat=60.0,
                    number=10,
                    radius=1000,
                    time=rn.start_time)
    o.run(steps=2)

    return o.plot(fast=False, buffer=.2, land_zorder=0)[1]
