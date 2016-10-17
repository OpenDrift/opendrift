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
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import unittest
from datetime import datetime, timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.pelagicegg import PelagicEggDrift


class TestRun(unittest.TestCase):

    def test_reader_boundary(self):
        o = PelagicEggDrift(loglevel=0)

        reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
        reader_nordic = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/fou-hi/nordic4km-1h/Nordic-4km_SURF_1h_avg_00.nc')
        o.add_reader([reader_norkyst, reader_nordic])
        o.fallback_values['land_binary_mask'] = 0
        lon = 0.2; lat = 61.0; # Close to NorKyst boundary
        time = reader_nordic.start_time
        o.seed_elements(lon, lat, radius=5000, number=100, time=time, z=0)
        o.config['turbulentmixing']['timestep'] = 5
        o.run(steps=5, time_step=3600,
              time_step_output=3600)


if __name__ == '__main__':
    unittest.main()
