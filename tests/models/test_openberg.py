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

from datetime import datetime, timedelta
import numpy as np
from opendrift.models.openberg import OpenBerg
from opendrift.readers import reader_netCDF_CF_generic

"""Tests for OpenBerg module."""


def test_openberg_constant_forcing():

    # Northwards current, no wind, no Coriolis
    o = OpenBerg(loglevel=50)
    o.set_config('environment:constant:x_sea_water_velocity', 0)
    o.set_config('environment:constant:y_sea_water_velocity', 1)
    o.set_config('environment:constant:x_wind', 0)
    o.set_config('environment:constant:y_wind', 0)
    o.set_config('drift:coriolis', False)
    o.seed_elements(4, 60, time=datetime.now())
    o.run(steps=2)
    np.testing.assert_almost_equal(o.result.lon[0][1], 4.0, 3)
    np.testing.assert_almost_equal(o.result.lat[0][1], 60.031, 3)

    # Northwards current, no wind, with Coriolis
    o = OpenBerg(loglevel=50)
    o.set_config('environment:constant:x_sea_water_velocity', 0)
    o.set_config('environment:constant:y_sea_water_velocity', 1)
    o.set_config('environment:constant:x_wind', 0)
    o.set_config('environment:constant:y_wind', 0)
    o.set_config('drift:coriolis', True)
    o.seed_elements(4, 60, time=datetime.now())
    o.run(steps=2)
    np.testing.assert_almost_equal(o.result.lon[0][1], 4.017, 3)
    np.testing.assert_almost_equal(o.result.lat[0][1], 60.031, 3)

    # Northwards current, eastwards wind, with Coriolis
    o = OpenBerg(loglevel=0)
    o.set_config('environment:constant:x_sea_water_velocity', 0)
    o.set_config('environment:constant:y_sea_water_velocity', 1)
    o.set_config('environment:constant:x_wind', 10)
    o.set_config('environment:constant:y_wind', 0)
    o.set_config('drift:coriolis', True)
    o.seed_elements(4, 60, time=datetime.now())
    o.run(steps=2)
    # The below is not expected: when adding eastwards wind, we expect iceberg to move
    # further eastwards, but we see that this is not the case
    np.testing.assert_almost_equal(o.result.lon[0][1], 4.0059, 3)
    np.testing.assert_almost_equal(o.result.lat[0][1], 60.0056, 3)


def test_openberg_norkyst():
    o = OpenBerg()
    reader_current = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
            '14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
    o.add_reader(reader_current)
    o.set_config('environment:constant:x_wind', 0)
    o.set_config('environment:constant:y_wind', 0)
    o.seed_elements(4., 62., time=reader_current.start_time)
    o.run(steps=3)

    np.testing.assert_almost_equal(o.result.lon[0][1], 3.998, 3)
    np.testing.assert_almost_equal(o.result.lat[0][1], 62.013, 3)
