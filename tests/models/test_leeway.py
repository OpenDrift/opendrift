#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import os
import time
from datetime import datetime, timedelta
from . import *

from opendrift.readers import reader_global_landmask
from opendrift.models.leeway import Leeway

"""Tests for Leeway module."""
def test_leewayprop():
    """Check that Leeway properties are properly read."""
    object_type = 85  # MED-WASTE-7
    lee = Leeway(loglevel=20)
    object_type = object_type
    assert lee.leewayprop[object_type]['Description'] == '>>Medical waste, syringes, small'
    assert lee.leewayprop[object_type]['DWSLOPE'] == 1.79

def test_leeway_config_object():
    """Check that correct object type is fetched from config"""
    l = Leeway(loglevel=20)
    l.set_config('seed:object_type', 'Surf board with person')
    l.set_config('environment:constant:x_wind', 0)
    l.set_config('environment:constant:y_wind', 0)
    l.set_config('environment:constant:x_sea_water_velocity', 0)
    l.set_config('environment:constant:y_sea_water_velocity', 0)
    l.seed_elements(lon=4.5, lat=60, number=100,
                    time=datetime(2015, 1, 1))
    objType = l.elements_scheduled.object_type
    assert l.leewayprop[objType]['Description'] == 'Surf board with person'
    assert l.leewayprop[objType]['OBJKEY'] == 'PERSON-POWERED-VESSEL-2'

def test_leewayrun(tmpdir, test_data):
    """Test the expected Leeway left/right split."""
    lee = Leeway(loglevel=20)
    object_type = 50  # FISHING-VESSEL-1
    reader_landmask = reader_global_landmask.Reader()
    lee.add_reader([reader_landmask])
    lee.set_config('general:coastline_approximation_precision', None)
    lee.set_config('environment:fallback:x_wind', 0)
    lee.set_config('environment:fallback:y_wind', 10)
    lee.set_config('environment:fallback:x_sea_water_velocity', 0)
    lee.set_config('environment:fallback:y_sea_water_velocity', 0)
    lee.seed_cone(lon=[4.5, 4.7], lat=[60.1, 60], number=100,
                  object_type=object_type,
                  time=[datetime(2015, 1, 1, 0), datetime(2015, 1, 1, 6)])
    # Check that 10 out of 100 elements strand towards coast
    lee.run(steps=24, time_step=3600)
    assert lee.num_elements_scheduled() == 0
    assert lee.num_elements_active() == 88
    assert lee.num_elements_deactivated() == 12  # stranded

    asciif = tmpdir + '/leeway_ascii.txt'
    lee.export_ascii(asciif)
    asciitarget = test_data + "/generated/test_leewayrun_export_ascii.txt"
    asciitarget2 = test_data + "/generated/test_leewayrun_export_ascii_v2.txt"
    from difflib import Differ
    with open(asciif) as file_1, open(asciitarget) as file_2:
        differ = Differ()
        for line in differ.compare(file_1.readlines(), file_2.readlines()):
            print(line)
    import filecmp
    if not filecmp.cmp(asciif, asciitarget):
        # Comparing with second versio of ASCII file, with slight numerical differences
        assert filecmp.cmp(asciif, asciitarget2)

def test_capsize():
    o = Leeway(loglevel=20)
    o.set_config('environment:constant', {'x_sea_water_velocity': 0, 'y_sea_water_velocity': 0,
                    'x_wind': 25, 'y_wind': 0, 'land_binary_mask': 0})
    o.set_config('processes:capsizing', True)
    o.set_config('capsizing:wind_threshold', 30)
    o.set_config('capsizing:wind_threshold_sigma', 3)
    o.set_config('capsizing:leeway_fraction', .4)
    o.seed_elements(lon=0, lat=60, time=datetime.now(), number=100)
    o.run(time_step=900, time_step_output=900, duration=timedelta(hours=6))
    assert o.elements.capsized.max() == 1
    assert o.elements.capsized.min() == 0
    assert o.elements.capsized.sum() == 18

    # Backward run, checking that forward capsizing is not happening
    ob = Leeway(loglevel=20)
    ob.set_config('processes:capsizing', True)
    ob.set_config('capsizing:wind_threshold', 30)
    ob.set_config('capsizing:wind_threshold_sigma', 3)
    ob.set_config('capsizing:leeway_fraction', .4)
    ob.set_config('environment:constant', {'x_sea_water_velocity': 0, 'y_sea_water_velocity': 0,
                    'x_wind': 25, 'y_wind': 0, 'land_binary_mask': 0})
    ob.seed_elements(lon=0, lat=60, time=datetime.now(), number=100)
    ob.run(time_step=-900, time_step_output=900, duration=timedelta(hours=6))
    assert ob.elements.capsized.max() == 0
    assert ob.elements.capsized.min() == 0
    assert ob.elements.capsized.sum() == 0

    # Backward run, checking that backward capsizing does happen
    ob = Leeway(loglevel=20)
    ob.set_config('processes:capsizing', True)
    ob.set_config('capsizing:wind_threshold', 30)
    ob.set_config('capsizing:wind_threshold_sigma', 3)
    ob.set_config('capsizing:leeway_fraction', .4)
    ob.set_config('environment:constant', {'x_sea_water_velocity': 0, 'y_sea_water_velocity': 0,
                    'x_wind': 25, 'y_wind': 0, 'land_binary_mask': 0})
    ob.seed_elements(lon=0, lat=60, time=datetime.now(), number=100, capsized=1)
    ob.run(time_step=-900, time_step_output=900, duration=timedelta(hours=6))
    assert ob.elements.capsized.max() == 1
    assert ob.elements.capsized.min() == 0
    assert ob.elements.capsized.sum() == 82
