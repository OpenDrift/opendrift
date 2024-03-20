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

import unittest
import pytest
from datetime import datetime, timedelta
import os
import inspect

import numpy as np

from opendrift.readers import reader_ArtificialOceanEddy
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_oscillating
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.basemodel import Mode, WrongMode
from opendrift.models.openoil import OpenOil
from opendrift.models.leeway import Leeway
from opendrift.models.pelagicegg import PelagicEggDrift
from opendrift.models.plastdrift import PlastDrift


class TestIO(unittest.TestCase):
    """Tests for file io"""

    def test_io_parquet(self):
        outfile = "test_io_parquet.nc"
        o = OceanDrift(
            loglevel=30,
            iomodule="parquet",
        )
        norkyst = reader_netCDF_CF_generic.Reader(
            o.test_data_folder()
            + "16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc"
        )
        landmask = reader_global_landmask.Reader()
        o.add_reader([landmask, norkyst])
        o.seed_elements(4.96, 60.1, radius=10, number=10, time=norkyst.start_time)
        o.run(
            steps=10,
            time_step=timedelta(minutes=30),
            time_step_output=timedelta(minutes=30),
            outfile=outfile,
            export_buffer_length=2,
        )
        os.remove(outfile)


if __name__ == "__main__":
    unittest.main()
