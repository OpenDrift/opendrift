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
import pytest
from datetime import datetime, timedelta
import numpy as np
import xarray as xr

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
try:
    import fastparquet
    has_fastparquet = True
except:
    has_fastparquet = False

need_fastparquet = pytest.mark.skipif(has_fastparquet == False,
        reason = 'fastparquet must be installed to use fastparquet writer')

def test_custom_result(tmpdir):
    """Adding custom data to self.result during self.prepare_run()"""

    def prepare_run(self):
        # Add custom scalar as attribute
        self.result['custom_scalar'] = 5
        # Add custom list as attribute
        self.result['custom_list'] = ['item1', 'item2', 'item3']
        # Add custom array as DataArray with dimensions other than (time, trajectory)
        da = xr.DataArray(
            data=np.random.rand(3, 3),
            dims=["dim1", "dim2"]
             )
        self.result['custom_array'] = da

    outfile = tmpdir + "test_custom_result.nc"
    o = OceanDrift(loglevel=50)
    time = datetime.now()
    o.prepare_run = prepare_run.__get__(o)  # Bind custom prepare_run to the instance
    o.seed_elements(lon=4, lat=60, time=time)
    o.run(steps=4, export_buffer_length=2, outfile=outfile)

    assert o.result.custom_array.shape == (3, 3)
    assert o.result.custom_scalar == 5
    assert len(o.result.custom_list) == 3
    ds = xr.open_dataset(outfile)
    assert(ds.custom_scalar == 5)
    assert 'time_coverage_end' in o.result.attrs
    assert 'time_coverage_end' in ds.attrs

@need_fastparquet
def test_io_parquet(tmpdir):
    outfile = tmpdir + "test_io_parquet.nc"
    o = OceanDrift(
        loglevel=30,
        iomodule="parquet",
    )
    norkyst = reader_netCDF_CF_generic.Reader(
        o.test_data_folder()
        + "16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc"
    )
    o.add_reader(norkyst)
    o.seed_elements(4.96, 60.1, radius=10, number=10, time=norkyst.start_time)
    o.run(
        steps=10,
        time_step=timedelta(minutes=30),
        time_step_output=timedelta(minutes=30),
        outfile=outfile,
        export_buffer_length=2,
    )
    os.remove(outfile)
