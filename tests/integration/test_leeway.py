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
from datetime import datetime

from opendrift.readers import reader_global_landmask
from opendrift.models.leeway import Leeway

"""Tests for Leeway module."""
def test_leewayprop():
    """Check that Leeway properties are properly read."""
    objectType = 85  # MED-WASTE-7
    lee = Leeway(loglevel=20)
    objectType = objectType
    assert lee.leewayprop[objectType]['Description'] == '>>Medical waste, syringes, small'
    assert lee.leewayprop[objectType]['DWSLOPE'] == 1.79

def test_leeway_config_object():
    """Check that correct object type is fetched from config"""
    l = Leeway(loglevel=20)
    l.set_config('seed:object_type', 'Surf board with person')
    l.seed_elements(lon=4.5, lat=60, number=100,
                    time=datetime(2015, 1, 1))
    objType = l.elements_scheduled.objectType
    assert l.leewayprop[objType]['Description'] == 'Surf board with person'
    assert l.leewayprop[objType]['OBJKEY'] == 'PERSON-POWERED-VESSEL-2'

def test_leewayrun(tmpdir):
    """Test the expected Leeway left/right split."""
    lee = Leeway(loglevel=30)
    objectType = 50  # FISHING-VESSEL-1
    reader_landmask = reader_global_landmask.Reader(extent=[ 3, 59.8, 6, 60.5 ])
    lee.add_reader([reader_landmask])
    lee.seed_elements(lon=4.5, lat=60, number=100,
                            objectType=objectType,
                            time=datetime(2015, 1, 1))
    lee.fallback_values['x_wind'] = 0
    lee.fallback_values['y_wind'] = 10
    lee.fallback_values['x_sea_water_velocity'] = 0
    lee.fallback_values['y_sea_water_velocity'] = 0
    # Check that 7 out of 100 elements strand towards coast
    lee.run(steps=24, time_step=3600)
    assert lee.num_elements_scheduled() == 0
    assert lee.num_elements_active() == 96
    assert lee.num_elements_deactivated() == 4  # stranded
    lee.export_ascii(tmpdir + '/leeway_ascii.txt')

