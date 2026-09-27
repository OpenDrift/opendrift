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
import numpy as np
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
    print('Comparing with first version of ASCII file')
    import filecmp
    if not filecmp.cmp(asciif, asciitarget):
        from difflib import Differ
        with open(asciif) as file_1, open(asciitarget) as file_2:
            differ = Differ()
            for line in differ.compare(file_1.readlines(), file_2.readlines()):
                print(line)
            # Comparing with second version of ASCII file, with slight numerical differences
            print('Comparing with second version of ASCII file')
            if not filecmp.cmp(asciif, asciitarget2):
                with open(asciif) as file_1, open(asciitarget2) as file_2:
                    differ = Differ()
                    for line in differ.compare(file_1.readlines(), file_2.readlines()):
                        print(line)
                raise ValueError('Leeway ascii output does not match any of the two template files')

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



def _run_constant_northward_wind(object_type, jibe_probability, hours=24):
    o = Leeway(loglevel=50)
    o.set_config('environment:constant', {'x_sea_water_velocity': 0, 'y_sea_water_velocity': 0,
                    'x_wind': 0, 'y_wind': 10, 'land_binary_mask': 0})
    o.seed_elements(lon=4, lat=60, time=datetime(2020, 1, 1), number=1000,
                    object_type=object_type, jibe_probability=jibe_probability)
    o.run(duration=timedelta(hours=hours), time_step=900)
    x_km = (o.elements.lon - 4) * 111.2 * np.cos(np.radians(60))
    y_km = (o.elements.lat - 60) * 111.2
    return o, x_km, y_km


def test_crosswind_direction(tmp_path):
    """Right-of-downwind elements should drift to the right of the wind.

    PIW-4 (survival suit) is asymmetric: right slope 1.36, offset -3.30,
    left slope -0.13, offset -2.65. With 10 m/s wind blowing towards north,
    right-drifting elements move east (~10.3 cm/s) and left-drifting west (~3.95 cm/s).
    A diagnostic plot is saved to tmp_path.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    o, x_km, y_km = _run_constant_northward_wind(object_type=4, jibe_probability=0)
    ori = o.elements.orientation
    p = o.leewayprop[4]
    to_km = .01 * 24 * 3600 / 1000  # cm/s over 24 hours
    expected = {0: (p['CWRSLOPE'] * 10 + p['CWROFFSET']) * to_km,
                1: (p['CWLSLOPE'] * 10 + p['CWLOFFSET']) * to_km}

    fig, ax = plt.subplots(figsize=(7, 7))
    for orientation, name, color in ((0, 'right', 'tab:red'), (1, 'left', 'tab:blue')):
        ind = ori == orientation
        ax.scatter(x_km[ind], y_km[ind], s=2, c=color, label=f'{name} of downwind (modelled)')
        ax.axvline(expected[orientation], color=color, ls='--',
                   label=f'{name}: expected mean {expected[orientation]:.1f} km')
        ax.plot(x_km[ind].mean(), y_km[ind].mean(), 'X', c=color, mec='k', ms=12)
    ax.plot(0, 0, 'k*', ms=10, label='seed')
    ax.annotate('wind', xy=(0, 5), xytext=(0, 0), ha='center',
                arrowprops=dict(arrowstyle='->', lw=2))
    ax.axvline(0, color='gray', lw=.5)
    ax.set_xlabel('Crosswind, east [km]')
    ax.set_ylabel('Downwind, north [km]')
    ax.set_title('PIW-4, 10 m/s wind towards north, 24 h, no jibing\n'
                 'X = modelled mean, dashed = expected from OBJECTPROP.DAT')
    ax.set_aspect('equal')
    leg = ax.legend(loc='lower left', fontsize=8)
    for h in leg.legend_handles:
        if hasattr(h, 'set_sizes'):
            h.set_sizes([20])
    plotfile = tmp_path / 'leeway_crosswind_direction.png'
    fig.savefig(plotfile, dpi=100)
    plt.close(fig)
    print(f'Plot saved to {plotfile}')

    assert y_km.mean() > 10  # Downwind is north
    np.testing.assert_allclose(x_km[ori == 0].mean(), expected[0], rtol=.1)  # right: east
    np.testing.assert_allclose(x_km[ori == 1].mean(), expected[1], rtol=.1)  # left: west


def test_jibing_swaps_crosswind_coefficients(tmp_path):
    """After jibing, elements must carry the coefficients of their new side.

    A diagnostic plot is saved to tmp_path.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    o, x_km, y_km = _run_constant_northward_wind(object_type=4, jibe_probability=.5, hours=6)
    p = o.leewayprop[4]
    ori = o.elements.orientation

    fig, ax = plt.subplots(figsize=(7, 6))
    for side, orientation, name, color in (('CWR', 0, 'right', 'tab:red'),
                                          ('CWL', 1, 'left', 'tab:blue')):
        ind = ori == orientation
        ax.scatter(o.elements.crosswind_slope[ind], o.elements.crosswind_offset[ind],
                   s=40, c=color, zorder=3, label=f'{ind.sum()} elements now {name} of downwind')
        ax.plot(p[side + 'SLOPE'], p[side + 'OFFSET'], 'o', mfc='none', mec=color, mew=2,
                ms=20, label=f'{side} coefficients in OBJECTPROP.DAT')
    ax.set_xlabel('Crosswind slope [%]')
    ax.set_ylabel('Crosswind offset [cm/s]')
    ax.set_title('PIW-4 after 6 h with jibe probability 0.5/h\n'
                 'each element should sit inside the circle of its side')
    ax.grid(True)
    ax.legend(fontsize=8)
    plotfile = tmp_path / 'leeway_jibing_coefficients.png'
    fig.savefig(plotfile, dpi=100)
    plt.close(fig)
    print(f'Plot saved to {plotfile}')

    assert np.any(ori == 0) and np.any(ori == 1)
    for side, orientation in (('CWR', 0), ('CWL', 1)):
        ind = ori == orientation
        np.testing.assert_allclose(o.elements.crosswind_slope[ind], p[side + 'SLOPE'])
        np.testing.assert_allclose(o.elements.crosswind_offset[ind], p[side + 'OFFSET'])
