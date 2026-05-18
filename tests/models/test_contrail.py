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

"""Tests for the ContrailDrift model."""

import pytest
import numpy as np
from datetime import datetime, timedelta

from opendrift.models.contrail import ContrailDrift, ContrailElement


@pytest.fixture
def contrail_dry():
    """ContrailDrift instance in completely dry air (immediate sublimation)."""
    c = ContrailDrift(loglevel=50)
    c.set_config('environment:constant:x_wind', 10)
    c.set_config('environment:constant:y_wind', 0)
    c.set_config('environment:constant:air_temperature', 220)
    c.set_config('environment:constant:specific_humidity', 0)   # completely dry
    c.set_config('environment:constant:horizontal_diffusivity', 0)
    return c


@pytest.fixture
def contrail_saturated():
    """ContrailDrift instance in saturated air (no sublimation)."""
    c = ContrailDrift(loglevel=50)
    c.set_config('environment:constant:x_wind', 10)
    c.set_config('environment:constant:y_wind', 0)
    c.set_config('environment:constant:air_temperature', 220)
    # Specific humidity that gives RHi >= 1 at 220 K, 11 km altitude
    # e_sat_ice(220) ≈ 3.5 Pa; p(11 km) ≈ 22 700 Pa
    # q_sat = 0.622 * e_sat / (p - e_sat) ≈ 0.622 * 3.5 / 22700 ≈ 9.6e-5 kg/kg
    c.set_config('environment:constant:specific_humidity', 1e-4)
    c.set_config('environment:constant:horizontal_diffusivity', 0)
    return c


# ------------------------------------------------------------------ #
# Basic instantiation and configuration                               #
# ------------------------------------------------------------------ #

def test_import():
    from opendrift.models.contrail import ContrailDrift
    assert ContrailDrift is not None


def test_contrail_element_variables():
    expected = {'ice_water_path', 'optical_depth', 'fall_speed'}
    assert expected.issubset(ContrailElement.variables.keys())


def test_coastline_action_default():
    c = ContrailDrift(loglevel=50)
    assert c.get_config('general:coastline_action') == 'none'


def test_ocean_only_default():
    c = ContrailDrift(loglevel=50)
    assert c.get_config('seed:ocean_only') is False


# ------------------------------------------------------------------ #
# Seeding                                                             #
# ------------------------------------------------------------------ #

def test_seed_above_land(contrail_dry):
    """Contrails can be seeded above land (seed:ocean_only = False)."""
    contrail_dry.seed_elements(
        lon=10.0, lat=48.0, z=11000,
        time=datetime(2024, 1, 1, 12))
    assert contrail_dry.num_elements_scheduled() == 1


def test_seed_multiple_elements(contrail_dry):
    lons = np.linspace(-10, 10, 50)
    lats = np.linspace(51, 53, 50)
    contrail_dry.seed_elements(
        lon=lons, lat=lats, z=11000,
        time=datetime(2024, 1, 1, 12))
    assert contrail_dry.num_elements_scheduled() == 50


def test_seed_custom_properties(contrail_saturated):
    contrail_saturated.seed_elements(
        lon=0., lat=50., z=10000,
        time=datetime(2024, 6, 1),
        ice_water_path=5e-4,
        optical_depth=0.5,
        fall_speed=0.05)
    contrail_saturated.run(steps=1, time_step=60)
    np.testing.assert_allclose(
        contrail_saturated.elements.ice_water_path, 5e-4, rtol=0.05)


# ------------------------------------------------------------------ #
# Advection                                                           #
# ------------------------------------------------------------------ #

def test_advection_with_wind(contrail_saturated):
    """Elements should move eastward with a 10 m/s zonal wind."""
    lon0, lat0 = 0.0, 50.0
    contrail_saturated.seed_elements(
        lon=lon0, lat=lat0, z=11000,
        time=datetime(2024, 1, 1))
    contrail_saturated.run(duration=timedelta(hours=1), time_step=300)

    lon_final = contrail_saturated.elements.lon[0]
    assert lon_final > lon0, 'Element should have moved eastward'


def test_no_advection_without_wind():
    """With zero wind, elements should not move horizontally."""
    c = ContrailDrift(loglevel=50)
    c.set_config('environment:constant:x_wind', 0)
    c.set_config('environment:constant:y_wind', 0)
    c.set_config('environment:constant:air_temperature', 220)
    c.set_config('environment:constant:specific_humidity', 1e-4)
    c.set_config('environment:constant:horizontal_diffusivity', 0)

    lon0, lat0 = 5.0, 55.0
    c.seed_elements(lon=lon0, lat=lat0, z=11000, time=datetime(2024, 1, 1))
    c.run(duration=timedelta(hours=2), time_step=300)

    np.testing.assert_allclose(c.elements.lon[0], lon0, atol=1e-4)
    np.testing.assert_allclose(c.elements.lat[0], lat0, atol=1e-4)


# ------------------------------------------------------------------ #
# Sedimentation                                                       #
# ------------------------------------------------------------------ #

def test_sedimentation_lowers_altitude(contrail_saturated):
    """Ice crystals should fall (z should decrease over time)."""
    z0 = 11000.
    contrail_saturated.seed_elements(
        lon=0., lat=50., z=z0, time=datetime(2024, 1, 1))
    contrail_saturated.run(duration=timedelta(hours=2), time_step=300)

    z_final = contrail_saturated.elements.z[0]
    assert z_final < z0, 'Element altitude should decrease due to sedimentation'


def test_sedimentation_rate(contrail_saturated):
    """Altitude decrease should match fall_speed × elapsed time."""
    z0 = 11000.
    fall_speed = 0.03   # m/s (default)
    duration_s = 3600.  # 1 hour

    contrail_saturated.seed_elements(
        lon=0., lat=50., z=z0, time=datetime(2024, 1, 1))
    contrail_saturated.run(duration=timedelta(seconds=duration_s),
                           time_step=300)

    expected_z = z0 - fall_speed * duration_s
    np.testing.assert_allclose(
        contrail_saturated.elements.z[0], expected_z, rtol=1e-3)


# ------------------------------------------------------------------ #
# Sublimation                                                         #
# ------------------------------------------------------------------ #

def test_sublimation_in_dry_air(contrail_dry):
    """Ice water path should decrease in completely dry air."""
    iwp0 = 1e-4
    contrail_dry.seed_elements(
        lon=0., lat=50., z=11000,
        time=datetime(2024, 1, 1),
        ice_water_path=iwp0)
    contrail_dry.run(duration=timedelta(minutes=30), time_step=60)

    iwp_final = contrail_dry.elements.ice_water_path
    # Allow for deactivation (empty array) or reduced IWP
    if len(iwp_final) > 0:
        assert np.all(iwp_final < iwp0), 'IWP should decrease in dry air'


def test_no_sublimation_in_saturated_air(contrail_saturated):
    """Ice water path should remain constant in fully saturated air."""
    iwp0 = 1e-4
    contrail_saturated.seed_elements(
        lon=0., lat=50., z=11000,
        time=datetime(2024, 1, 1),
        ice_water_path=iwp0)
    contrail_saturated.run(duration=timedelta(hours=2), time_step=300)

    np.testing.assert_allclose(
        contrail_saturated.elements.ice_water_path[0], iwp0, rtol=1e-4)


def test_elements_deactivated_after_sublimation(contrail_dry):
    """All elements should be deactivated after sufficient time in dry air."""
    contrail_dry.set_config('contrail:sublimation_timescale', 300.)
    contrail_dry.seed_elements(
        lon=np.linspace(0, 5, 10),
        lat=np.full(10, 50.),
        z=11000,
        time=datetime(2024, 1, 1),
        ice_water_path=1e-4)
    contrail_dry.run(duration=timedelta(hours=3), time_step=60)

    assert contrail_dry.num_elements_active() == 0, (
        'All elements should sublimate in completely dry air')
    assert contrail_dry.num_elements_total() == 10


# ------------------------------------------------------------------ #
# Thermodynamics helper
# ------------------------------------------------------------------ #

def test_saturation_pressure_over_ice():
    """e_sat_ice should be ~611 Pa at 273.16 K (triple point)."""
    e = ContrailDrift._saturation_pressure_over_ice(np.array([273.16]))
    np.testing.assert_allclose(e, 611.73, rtol=1e-3)


def test_saturation_pressure_increases_with_temperature():
    T = np.array([200., 220., 240., 260.])
    e = ContrailDrift._saturation_pressure_over_ice(T)
    assert np.all(np.diff(e) > 0), 'e_sat should increase with temperature'
