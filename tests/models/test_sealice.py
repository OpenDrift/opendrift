import numpy as np
import pytest
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.sealice import SeaLice
from datetime import timedelta


@pytest.mark.xfail(strict=False)
def test_sealice_larc():
    o = SeaLice(loglevel=30)
    reader_arome = reader_netCDF_CF_generic.Reader(
        o.test_data_folder() +
        '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF_20160202_subset.nc')
    o.add_reader([reader_arome])
    reader_light = reader_netCDF_CF_generic.Reader(
        'https://opendap.larc.nasa.gov/opendap/hyrax/POWER/monthly/power_901_monthly_radiation_utc.nc',
        standard_name_mapping={
            'ALLSKY_SFC_SW_DWN': 'surface_net_downward_radiative_flux'
        })
    o.add_reader(reader_light)
    lcts = timedelta(hours=1).total_seconds()  #seeding time-steps
    lat = 67.711251
    lon = 13.556971  # Lofoten
    o.set_config('lice:seeding_time_step', lcts)
    o.seed_elements(lon,
                    lat,
                    radius=5000,
                    number=1000,
                    time=reader_arome.start_time,
                    particle_biomass=2000,
                    z=-5)
    o.run(time_step=lcts, steps=2)
    np.testing.assert_almost_equal(o.elements.lon.max(), 13.882, 2)
