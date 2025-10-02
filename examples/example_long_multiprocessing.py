#!/usr/bin/env python
"""
Multiprocessing - Parallel execution of multiple OpenDrift instances in distinct CPUs
======================================================================================
"""

from datetime import datetime
import numpy as np
import xarray as xr
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
from multiprocessing import Pool
import glob
import opendrift

def concatenate_outputs(input_files, output_file):
    """
    Concatenate OpenDrift NetCDF outputs from multiple workers into one single output file.

    """

    ds_list = []
    trajectory_offset = 0
    for f in input_files:
            dsi = xr.open_dataset(f)
            traj_vars = [var for var in list(dsi.variables) if "trajectory" in dsi[var].dims]
            ds_traj = dsi[traj_vars]
            ds_traj = ds_traj.assign_coords(trajectory=ds_traj.trajectory + trajectory_offset)

            trajectory_offset += ds_traj.dims["trajectory"]
            ds_list.append(ds_traj)

    ds_conc = xr.concat(ds_list, dim="trajectory")
    ds_conc.to_netcdf(output_file)

    return ds_conc

def RunOceanDrift(pool_number):
    """
    Configure and run one OpenDrift instance

    """

    o = OceanDrift(loglevel=0, logfile="output_W"+str(pool_number)+".log", seed=0)

    reader_topaz4 = reader_netCDF_CF_generic.Reader("https://thredds.met.no/thredds/dodsC/topaz/dataset-topaz4-arc-myoceanv2-be")
    o.add_reader(reader_topaz4, variables=['x_sea_water_velocity', 'y_sea_water_velocity','sea_water_temperature','sea_water_salinity','sea_floor_depth_below_sea_level'])

    o.set_config('drift:horizontal_diffusivity', 50)

    time = datetime(2023,1,1)

    ntraj=500
    iniz=np.random.rand(ntraj) * -10. # seeding the chemicals in the upper 10m

    o.seed_elements(lat=positions[pool_number][0], lon=positions[pool_number][1], z=iniz, radius=2000, number=ntraj, time=time, origin_marker=np.ones(ntraj)*(pool_number))

    o.run(steps=7*4, time_step=3600*6, time_step_output=3600*6, outfile = "output_W"+str(pool_number)+".nc")


positions=[(58.5,3),(58.2, 2.5),(58,1),(57.8,1.2),(57,1),(56.8,1.8),(56.5,2),(56,2.2)]
pool_size=len(positions)

#%% Run pool of OpenDrift instances in parallel using distinct CPUs and concatenate results
with Pool(pool_size) as p:
    p.starmap(RunOceanDrift, [(i,) for i in range(pool_size)])

concatenate_outputs(input_files=glob.glob("output_W*.nc"), output_file="output_total.nc")


#%% Generates animation of the first instance and the concatenated output

o0 = opendrift.open("output_W0.nc")

#%%
# .. image:: /gallery/animations/example_multiprocessing_0.gif

o0.animation(color='origin_marker',
            markersize=3,
            vmin=0,vmax=pool_size-1,
            colorbar=False,
            fast = True,
            lscale = 'l')

o = opendrift.open("output_total.nc")

#%%
# .. image:: /gallery/animations/example_multiprocessing_1.gif

o.animation(color='origin_marker',
            markersize=3,
            vmin=0,vmax=pool_size-1,
            colorbar=False,
            fast = True,
            lscale = 'l')
