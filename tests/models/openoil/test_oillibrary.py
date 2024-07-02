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
import unittest
from datetime import datetime, timedelta

import xarray as xr
import numpy as np

from opendrift.models.openoil import OpenOil


def test_new_oil():
    o = OpenOil(loglevel=50, location='Norway')
    oiltypes = o._config['seed:oil_type']['enum']
    assert 'HEIDRUN AARE 2023' in oiltypes
    assert len(oiltypes) >= 179

def test_oils():
    o = OpenOil(loglevel=50, weathering_model='noaa')

    assert len(o.oiltypes) >= 1247

    for oiltype in o.oiltypes[12:14]:
        if oiltype == 'JP-8':
            continue
        o = OpenOil(loglevel=50, weathering_model='noaa')
        o.set_config('environment:fallback:x_wind', 7)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', .7)
        o.set_config('environment:fallback:y_sea_water_velocity', 0)
        o.set_config('processes:evaporation', True)
        o.set_config('processes:emulsification', True)
        o.set_config('drift:vertical_mixing', False)
        o.set_config('drift:wind_uncertainty', 0)
        o.set_config('drift:current_uncertainty', 0)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.seed_elements(4.7,
                        60.0,
                        radius=3000,
                        number=3,
                        z=0,
                        time=datetime.now(),
                        oil_type=oiltype)
        o.run(steps=3)
        initial_mass = o.get_property('mass_oil')[0][0, 0]
        assert o.elements.mass_evaporated.min(
        ) == o.elements.mass_evaporated.max()
        assert o.elements.mass_evaporated.min() > 0
        assert o.elements.mass_evaporated.max() / initial_mass <= 1
        print(oiltype, o.elements.mass_evaporated.min())

def test_oilbudget():
    for windspeed in [3, 8]:
        for dispersion in [True, False]:
            o = OpenOil(loglevel=30, weathering_model='noaa')
            o.set_config('environment:fallback:x_wind', windspeed)
            o.set_config('environment:fallback:y_wind', 0)
            o.set_config('environment:fallback:x_sea_water_velocity', 0)
            o.set_config('environment:fallback:y_sea_water_velocity', 0)
            o.set_config('environment:fallback:land_binary_mask', 0)
            o.set_config('seed:oil_type', 'SIRTICA')
            o.set_config('processes:dispersion', dispersion)
            seed_hours = 3
            m3_per_hour = 200
            o.set_config('seed:m3_per_hour', m3_per_hour)
            o.seed_elements(
                lon=0,
                lat=60,
                number=100,
                #m3_per_hour=m3_per_hour,
                time=[
                    datetime.now(),
                    datetime.now() + timedelta(hours=seed_hours)
                ])
            o.run(duration=timedelta(hours=4))
            b = o.get_oil_budget()
            density = o.get_property('density')[0][0, 0]
            volume = b['mass_total'] / density
            np.testing.assert_almost_equal(volume[-1],
                                           seed_hours * m3_per_hour, 2)

            if dispersion is True:
                disp = 'dispersion'
            else:
                disp = 'nodispersion'
            filename = 'oilbudget_%s_%s.png' % (windspeed, disp)
            o.plot_oil_budget(filename=filename)
            os.remove(filename)

def test_oil_volume():
    o = OpenOil(loglevel=50)
    m3_per_hour = 50
    o.set_config('seed:m3_per_hour', m3_per_hour)
    o.set_config('environment:fallback:x_wind', 0)
    o.set_config('environment:fallback:y_wind', 0)
    o.set_config('environment:fallback:x_sea_water_velocity', 0)
    o.set_config('environment:fallback:y_sea_water_velocity', 0)
    o.set_config('environment:fallback:land_binary_mask', 0)
    o.seed_elements(lon=4, lat=60, time=datetime.now())
    o.run(steps=1)
    b = o.get_oil_budget()
    volume = b['mass_total'] / b['oil_density']
    assert volume[0] == m3_per_hour
    assert volume[1] == m3_per_hour


def test_dispersion():
    # Condensate is temporarily not gnome_suitable, cuts shuld be added
    #for oil in ['SMORBUKK KONDENSAT', 'SKRUGARD']:
    for oil in ['SKRUGARD']:
        for windspeed in [3, 8]:
            if oil == 'SKRUGARD' and windspeed == 3:
                continue
            o = OpenOil(loglevel=20, weathering_model='noaa')
            o.set_config('processes:dispersion', True)
            o.set_config('vertical_mixing:timestep', 10)
            o.set_config('environment:fallback:land_binary_mask', 0)
            o.set_config('environment:fallback:x_wind', windspeed)
            o.set_config('environment:fallback:y_wind', 0)
            o.set_config('environment:fallback:x_sea_water_velocity', 0)
            o.set_config('environment:fallback:y_sea_water_velocity', .3)
            oilname = oil
            if oil not in o.oiltypes:
                if oilname == 'SKRUGARD':
                    oilname = 'SKRUGARD 2012'
                elif oilname == 'SMORBUKK KONDENSAT':
                    oilname = 'SMORBUKK KONDENSAT 2003'
            o.seed_elements(lon=4.8,
                            lat=60,
                            number=100,
                            time=datetime.now(),
                            oil_type=oilname)
            o.run(duration=timedelta(hours=3), time_step=900)

            b = o.get_oil_budget()
            actual_dispersed = b['mass_dispersed'] / b['mass_total']
            actual_submerged = b['mass_submerged'] / b['mass_total']
            actual_evaporated = b['mass_evaporated'] / b['mass_total']
            print('Dispersion fraction %f for '
                  '%s and wind speed %f' %
                  (actual_dispersed[-1], oil, windspeed))
            if oil == 'SMORBUKK KONDENSAT' and windspeed == 3:
                fraction_dispersed = 0
                fraction_submerged = 0
                fraction_evaporated = 0.589
                meanlon = 4.81742
            elif oil == 'SMORBUKK KONDENSAT' and windspeed == 8:
                fraction_dispersed = 0.092
                fraction_submerged = 0.288
                fraction_evaporated = 0.517
                meanlon = 4.816
            elif oil == 'SKRUGARD' and windspeed == 8:
                fraction_dispersed = 0.139
                fraction_submerged = 0.418
                fraction_evaporated = 0.123
                meanlon = 4.825
            else:
                fraction_dispersed = -1  # not defined

            np.testing.assert_almost_equal(actual_dispersed[-1],
                                           fraction_dispersed, 2)
            np.testing.assert_almost_equal(actual_submerged[-1],
                                           fraction_submerged, 2)
            np.testing.assert_almost_equal(actual_evaporated[-1],
                                           fraction_evaporated, 2)
            np.testing.assert_almost_equal(np.mean(o.elements.lon), meanlon, 3)
            #o.plot_oil_budget()


def test_no_dispersion():
    o = OpenOil(loglevel=50, weathering_model='noaa')

    o.set_config('processes:dispersion', False)
    o.set_config('environment:fallback:land_binary_mask', 0)
    o.set_config('environment:fallback:x_wind', 8)
    o.set_config('environment:fallback:y_wind', 8)
    o.set_config('environment:fallback:x_sea_water_velocity', 0)
    o.set_config('environment:fallback:y_sea_water_velocity', .3)
    o.seed_elements(lon=4.8,
                    lat=60,
                    number=100,
                    time=datetime.now(),
                    oil_type='SIRTICA')
    o.run(duration=timedelta(hours=2), time_step=1800)
    b = o.get_oil_budget()
    actual_dispersed = b['mass_dispersed'] / b['mass_total']
    np.testing.assert_almost_equal(actual_dispersed[-1], 0)
    np.testing.assert_array_almost_equal(o.elements.lon[4:7], [4.816, 4.797, 4.803], 3)


def test_biodegradation():
    o = OpenOil(loglevel=50, weathering_model='noaa')

    o.set_config('processes:dispersion', True)
    o.set_config('processes:evaporation', True)
    o.set_config('processes:emulsification', True)
    o.set_config('processes:biodegradation', True)
    o.set_config('environment:fallback:land_binary_mask', 0)
    o.set_config('environment:fallback:x_wind', 0)
    o.set_config('environment:fallback:y_wind', 0)
    o.set_config('environment:fallback:x_sea_water_velocity', 0)
    o.set_config('environment:fallback:y_sea_water_velocity', 0)
    o.set_config('environment:fallback:sea_water_temperature', 30)
 
    o.seed_elements(lon=4.8,
                    lat=60,
                    number=100,
                    time=datetime.now(),
                    oil_type='SIRTICA')
    o.run(duration=timedelta(days=1), time_step=1800)
    initial_mass = o.get_property('mass_oil')[0][0, 0]
    biodegraded30 = o.elements.mass_biodegraded
    factor = 0.120  #(1-e^(-1))
    np.testing.assert_almost_equal(biodegraded30[-1] / initial_mass, factor, 3)


def test_droplet_distribution():
    # TODO: Should also add Li2017 here
    for droplet_distribution in ['Johansen et al. (2015)']:
        o = OpenOil(loglevel=50, weathering_model='noaa')
        if 'SKRUGARD' in o.oiltypes:
            oiltype = 'SKRUGARD'
        else:
            oiltype = 'SKRUGARD 2012'
        o.set_config('wave_entrainment:droplet_size_distribution', droplet_distribution)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('environment:fallback:x_wind', 8)
        o.set_config('environment:fallback:y_wind', 0)
        o.set_config('environment:fallback:x_sea_water_velocity', 0)
        o.set_config('environment:fallback:y_sea_water_velocity', .3)
        o.seed_elements(lon=4.8,
                        lat=60,
                        number=100,
                        time=datetime.now(),
                        oil_type=oiltype)
        o.run(duration=timedelta(hours=1), time_step=1800)
        d = o.elements.diameter
        # Suspicious, Sintef-param should give larer droplets
        if droplet_distribution == 'Johansen et al. (2015)':
            #self.assertAlmostEqual(d.mean(), 0.000072158)
            np.testing.assert_almost_equal(d.mean(), 0.000653, 2)


def test_seed_metadata_output():
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    outfile = 'openoil.nc'
    a = {
        'lon': 4.8,
        'lat': 60,
        'z': -10,
        'number': 10,
        'radius': 100,
        'time': datetime(2021, 1, 22)
    }
    o = OpenOil(loglevel=50)
    o.set_config('environment:fallback:land_binary_mask', 0)
    o.set_config('environment:fallback:x_wind', 1)
    o.set_config('environment:fallback:y_wind', 0)
    o.set_config('environment:fallback:x_sea_water_velocity', 0)
    o.set_config('environment:fallback:y_sea_water_velocity', .3)
    o.seed_elements(lon=a['lon'],
                    lat=a['lat'],
                    z=a['z'],
                    radius=a['radius'],
                    number=a['number'],
                    time=a['time'])
    o.run(duration=timedelta(hours=1), time_step=1800, outfile=outfile)
    f = xr.open_dataset(outfile)
    for var in ['lon', 'lat', 'number', 'radius', 'z', 'm3_per_hour']:
        s = 'seed_' + var
        if var not in a:
            a[var] = o.get_config('seed:' + var)
        np.testing.assert_almost_equal(float(a[var]), float(f.attrs[s]), 3)
    assert str(a['time']) == f.attrs['seed_time']
    for var in [
            'seed_m3_per_hour', 'time_coverage_start', 'time_coverage_end',
            'seed_oiltype', 'simulation_time'
    ]:
        assert var in f.attrs
    f.close()
    os.remove(outfile)
    # Same for cone-seeding
    a = {
        'lon': [4.7, 4.8],
        'lat': [60, 60.1],
        'number': 10,
        'radius': [0, 100],
        'time': [datetime(2021, 1, 22, 0),
                 datetime(2021, 1, 22, 6)]
    }
    a_config = {'seafloor': True, 'm3_per_hour': 55}
    o = OpenOil(loglevel=0)
    o.set_config('environment:fallback:land_binary_mask', 0)
    o.set_config('environment:fallback:x_wind', 1)
    o.set_config('environment:fallback:y_wind', 0)
    o.set_config('environment:fallback:x_sea_water_velocity', 0)
    o.set_config('environment:fallback:y_sea_water_velocity', .3)
    o.set_config('environment:fallback:y_sea_water_velocity', .3)
    for ac in a_config:
        o.set_config('seed:' + ac, a_config[ac])
    o.seed_cone(lon=a['lon'],
                lat=a['lat'],
                radius=a['radius'],
                number=a['number'],
                time=a['time'])
    o.seed_cone(lon=a['lon'],
                lat=a['lat'],
                radius=a['radius'],
                number=a['number'] + 4,
                time=a['time'])
    o.run(duration=timedelta(hours=1), time_step=1800, outfile=outfile)
    f = xr.open_dataset(outfile)
    import geojson
    gj = geojson.loads(f.attrs['seed_geojson'])
    print(gj)
    fe = gj['features']
    for var in ['lon', 'lat', 'number', 'radius', 'm3_per_hour']:
        s = 'seed_' + var
        if var not in a:
            a[var] = o.get_config('seed:' + var)
        np.testing.assert_almost_equal(
            np.atleast_1d(a[var]).mean(), float(f.attrs[s]), 3)
        if var not in ['lon', 'lat']:  # Check also versus GeoJSON attribute
            geojsonvar = np.atleast_1d(fe[0]['properties'][var]).mean()
            avar = np.atleast_1d(a[var]).mean()
            assert avar == geojsonvar

    assert fe[0]['properties']['z'] == 'seafloor'
    assert str(a['time'][0]) == f.attrs['seed_time']
    for var in [
            'seed_m3_per_hour', 'time_coverage_start', 'time_coverage_end',
            'seed_oiltype', 'simulation_time'
    ]:
        assert var in f.attrs
    os.remove(outfile)
