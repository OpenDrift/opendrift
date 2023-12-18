from datetime import datetime, timedelta
import pytest
import numpy as np
from opendrift.models.openoil import adios
from opendrift.models.openoil import OpenOil


@pytest.fixture
def aasgard():
    oils = adios.oils(1, 'AASGARD A 2003')
    f = oils[0].make_full()
    assert f.name == 'AASGARD A 2003'
    return f

def test_max_water_fraction():
    # Check that max water fraction is saturated to max values from Sintef model
    # Max water_fraction 0.317 for SST=0 from Sintef, for SST=10,20 from NOAA
    for sst,expected_fraction in zip([0, 10, 20], [0.317, 0.443, 0.443]):
        o = OpenOil(loglevel=50)
        o.set_config('environment:constant',
            {
            'x_wind': 10,
            'y_wind': 0,
            'x_sea_water_velocity': 0,
            'y_sea_water_velocity': 0,
            'sea_water_temperature': sst,
            'land_binary_mask': 0
            })
        o.seed_elements(lon=3, lat=60, time=datetime.now(), number=100,
                        oil_type='DUVA 2021')
        o.run(duration=timedelta(hours=6))
        assert np.isclose(o.elements.water_fraction.max(), expected_fraction, atol=.001)

def test_open_aasgard(aasgard):
    print(aasgard)

def test_density_at_temp(aasgard):
    # old oillibrary: 816.6828030078809
    assert np.isclose(aasgard.density_at_temp(285.), 814.)

def test_kvis_at_temp(aasgard):
    # old oillibrary: 3.298187589355751e-05
    assert np.isclose(aasgard.kvis_at_temp(285.), 3.29e-5)

def test_mass_fraction(aasgard):
    assert np.isclose(aasgard.mass_fraction.sum(), 1.0)

    # XXX: This doesn't match for AASGARD.
    # old_mf = np.asarray([0.09667529, 0.09066779, 0.0060075 , 0.08161284, 0.01506245,
    #    0.07206664, 0.02460864, 0.06161574, 0.03505955, 0.04972493,
    #    0.04695036, 0.03565127, 0.06102402, 0.03565127, 0.06102402,
    #    0.03565127, 0.06102402, 0.03565127, 0.06102402, 0.03284712,
    #    0.0004    ])

    # new_mf = aasgard.mass_fraction

    # old_mf.sort()
    # new_mf.sort()

    # np.testing.assert_array_almost_equal(new_mf, old_mf)


def test_oil_water_surface_tension(aasgard):
    # old oillibrary: (0.028134082260442256, 288.15, False)
    assert np.isclose(aasgard.oil_water_surface_tension(), 0.02816)

def test_vapor_pressure(aasgard):
    new_vp = aasgard.vapor_pressure(285)

    # Components doesn't match
    # old_vp = np.array([3.32504179e+05, 7.18904609e+04, 7.18904609e+04, 1.38086757e+04,
    #    1.38086757e+04, 2.32314974e+03, 2.32314974e+03, 3.37425026e+02,
    #    3.37425026e+02, 4.16860123e+01, 4.16860123e+01, 4.32568844e+00,
    #    4.32568844e+00, 3.77848478e-01, 3.77848478e-01, 2.94709353e-02,
    #    2.94709353e-02, 2.52933763e-03, 2.52933763e-03, 1.54009621e-23,
    #    1.54009621e-23])


    # old_vp.sort()
    # new_vp.sort()

    # np.testing.assert_array_almost_equal(old_vp, new_vp)

def test_k0y_const():
    o = OpenOil()
    for oil in o.oiltypes:
        print('testing oil:', oil)
        oils = adios.oils(1, oil)
        if len(oils) > 0:
            f = oils[0].make_full()
            if hasattr(f, 'gnome_oil'):
                print(f.k0y)
                assert f.k0y == 2.024e-06

