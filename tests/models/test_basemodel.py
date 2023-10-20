from datetime import datetime, timedelta
from opendrift.readers import reader_global_landmask
from opendrift.models.leeway import Leeway
import numpy as np
import pytest

def test_simulation_back_extent():
    # backward
    leeb = Leeway()

    leeb.set_config('environment:fallback:x_wind', 1.5)
    leeb.set_config('environment:fallback:y_wind', 10)
    leeb.set_config('environment:fallback:x_sea_water_velocity', 1.5) # maximum speed in automatic landmask
    leeb.set_config('environment:fallback:y_sea_water_velocity', 0)

    object_type = 50  # FISHING-VESSEL-1
    leeb.seed_elements(lon=4, lat=60, number=100,
                            object_type=object_type,
                            time=datetime(2015, 1, 1))

    with pytest.raises(ValueError) as ex:
        leeb.run(duration=timedelta(days=-1), time_step=3600, time_step_output=10*3600)
    assert 'Time step must be negative if duration is negative.' in str(ex.value)

def test_simulation_matches_forw_backward():
    """ Check if simulation extent matches for both forward and backward modeling. """
    # forward
    leef = Leeway()

    leef.set_config('environment:fallback:x_wind', -1.5)
    leef.set_config('environment:fallback:y_wind', -10)
    leef.set_config('environment:fallback:x_sea_water_velocity', -1.5) # maximum speed in automatic landmask
    leef.set_config('environment:fallback:y_sea_water_velocity', 0)

    object_type = 50  # FISHING-VESSEL-1
    leef.seed_elements(lon=4.5, lat=60, number=100,
                            object_type=object_type,
                            time=datetime(2015, 1, 1))

    leef.run(steps=2, time_step=10*3600, time_step_output=10*3600)

    # backward
    leeb = Leeway()

    leeb.set_config('environment:fallback:x_wind', 1.5)
    leeb.set_config('environment:fallback:y_wind', 10)
    leeb.set_config('environment:fallback:x_sea_water_velocity', 1.5) # maximum speed in automatic landmask
    leeb.set_config('environment:fallback:y_sea_water_velocity', 0)

    object_type = 50  # FISHING-VESSEL-1
    leeb.seed_elements(lon=4.5, lat=60, number=100,
                            object_type=object_type,
                            time=datetime(2015, 1, 1))

    leeb.run(steps=2, time_step=-10*3600, time_step_output=-10*3600)

    maskf = leef.env.readers['global_landmask']
    maskb = leeb.env.readers['global_landmask']
    assert maskf.xmin == maskb.xmin
    assert maskf.ymin == maskb.ymin
    assert maskf.xmax == maskb.xmax
    assert maskf.ymax == maskb.ymax

    assert leef.num_elements_scheduled() == leeb.num_elements_scheduled()
    assert leef.num_elements_active() == leeb.num_elements_active()
    assert leef.num_elements_deactivated() == leeb.num_elements_deactivated()

    flon, flat = leef.get_lonlats()
    flon = flon[:,-1]
    flat = flat[:,-1]

    blon, blat = leeb.get_lonlats()
    blon = blon[:,-1]
    blat = blat[:,-1]

    np.testing.assert_array_almost_equal(np.sort(flon), np.sort(blon))
    np.testing.assert_array_almost_equal(np.sort(flat), np.sort(blat), decimal = 5)

def test_clone():
    lee = Leeway()
    lee.set_config('environment:fallback:x_wind', 1.5)
    lee.set_config('environment:fallback:y_wind', 10)
    lee.set_config('environment:fallback:x_sea_water_velocity', 1.5) # maximum speed in automatic landmask
    lee.set_config('environment:fallback:y_sea_water_velocity', 0)

    object_type = 50  # FISHING-VESSEL-1
    lee.seed_elements(lon=4.5, lat=60, number=100,
                            object_type=object_type,
                            time=datetime(2015, 1, 1))
    lee.run(steps=2, time_step=-10*3600, time_step_output=-10*3600)

    assert lee.get_config('environment:fallback:x_wind') == 1.5

    c = lee.clone()
    assert isinstance(c, Leeway)
    assert c.env.readers == lee.env.readers
    assert len(c.elements) == 0
    assert c.get_config('environment:fallback:x_wind') == 1.5

    object_type = 50  # FISHING-VESSEL-1
    c.seed_elements(lon=4.5, lat=60, number=100,
                            object_type=object_type,
                            time=datetime(2015, 1, 1))
    c.run(steps=2, time_step=-10*3600, time_step_output=-10*3600)

