import sys
import logging
from datetime import datetime, timedelta
from opendrift.readers import reader_global_landmask
from opendrift.models.leeway import Leeway
from opendrift.models.oceandrift import OceanDrift
import numpy as np
import pytest

def test_logging(tmpdir, capsys):
    # Accepting small variations in log output,
    # depending on machine, and from which folder test is run
    accepted = (288, 291, 316) 

    # Logging to console
    logfile = None
    o = OceanDrift(loglevel=0, logfile=logfile)
    o.set_config('environment:fallback:x_wind', 1.5)
    o.set_config('environment:fallback:y_wind', 10)
    o.set_config('environment:fallback:x_sea_water_velocity', .5)
    o.set_config('environment:fallback:y_sea_water_velocity', 0)
    o.seed_elements(lon=4, lat=60, number=5, time=datetime(2025, 1, 1))
    o.run(steps=3, time_step=600, time_step_output=600)
        
    stdout, stderr = capsys.readouterr()
    assert 'step 3 of 3 - 5 active elements (0 deactivated)' in stderr
    assert 'Changed mode from Mode.Run to Mode.Result' in stderr
    assert len(stderr.splitlines()) in accepted  # Depending on from where pytest is run

    # Logging to file
    logfile = tmpdir+'/test_log.txt'
    o = OceanDrift(loglevel=0, logfile=logfile)
    o.set_config('environment:fallback:x_wind', 1.5)
    o.set_config('environment:fallback:y_wind', 10)
    o.set_config('environment:fallback:x_sea_water_velocity', .5)
    o.set_config('environment:fallback:y_sea_water_velocity', 0)
    o.seed_elements(lon=4, lat=60, number=5, time=datetime(2025, 1, 1))
    o.run(steps=3, time_step=600, time_step_output=600)
    log = open(logfile).read()
    assert 'step 3 of 3 - 5 active elements (0 deactivated)' in log
    assert 'Changed mode from Mode.Run to Mode.Result' in log
    assert len(log.splitlines()) in accepted

    # Logging to existing file, assure overwrite and not append
    logfile = tmpdir+'/test_log.txt'
    o = OceanDrift(loglevel=0, logfile=logfile)
    o.set_config('environment:fallback:x_wind', 1.5)
    o.set_config('environment:fallback:y_wind', 10)
    o.set_config('environment:fallback:x_sea_water_velocity', .5)
    o.set_config('environment:fallback:y_sea_water_velocity', 0)
    o.seed_elements(lon=4, lat=60, number=5, time=datetime(2025, 1, 1))
    o.run(steps=3, time_step=600, time_step_output=600)
    log = open(logfile).read()
    assert 'step 3 of 3 - 5 active elements (0 deactivated)' in log
    assert 'Changed mode from Mode.Run to Mode.Result' in log
    assert len(log.splitlines()) in accepted

    # Logging to both file and terminal
    logfile = [tmpdir+'/test_log2.txt', logging.StreamHandler(sys.stdout)]
    o = OceanDrift(loglevel=0, logfile=logfile)
    o.set_config('environment:fallback:x_wind', 1.5)
    o.set_config('environment:fallback:y_wind', 10)
    o.set_config('environment:fallback:x_sea_water_velocity', .5)
    o.set_config('environment:fallback:y_sea_water_velocity', 0)
    o.seed_elements(lon=4, lat=60, number=5, time=datetime(2025, 1, 1))
    o.run(steps=3, time_step=600, time_step_output=600)
    log = open(logfile[0]).read()
    assert 'step 3 of 3 - 5 active elements (0 deactivated)' in log
    assert 'Changed mode from Mode.Run to Mode.Result' in log
    stdout, stderr = capsys.readouterr()
    assert 'step 3 of 3 - 5 active elements (0 deactivated)' in stdout
    assert 'Changed mode from Mode.Run to Mode.Result' in stdout
    assert len(stdout.splitlines()) in accepted
    assert len(log.splitlines()) in accepted
    assert stderr != log  # Since times are different


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

    flon = leef.result.lon.isel(time=-1)
    flat = leef.result.lat.isel(time=-1)
    blon = leeb.result.lon.isel(time=-1)
    blat = leeb.result.lat.isel(time=-1)

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

