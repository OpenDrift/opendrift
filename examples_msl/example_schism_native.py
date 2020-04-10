#!/usr/bin/env python
"""
SCHISM native reader
==================================
"""
# docker run -it  --rm -v C:/Users/Simon/Documents/GitHub/opendrift_simon:/code/ -v E:/metocean/0472_SpatNZ_MarlboroughSounds/schism_flows_netcdf/:/data opendrift/opendrift:py3-v1.0.7

import numpy as np
from opendrift.readers import reader_netCDF_CF_unstructured_SCHISM
from opendrift.readers import reader_global_landmask
from opendrift.models.oceandrift import OceanDrift

###############################
# MODEL
###############################
o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
###############################
# READERS
###############################
# Creating and adding reader using a native SCHISM netcdf output file
# SCHISM reader
# Landmask - (uses cartopy+rasterized GHSS shorelines)
reader_landmask = reader_global_landmask.Reader(
                    llcrnrlon=171.5, llcrnrlat=-43.5,
                    urcrnrlon=177.0, urcrnrlat=-38.0)
# NZTM proj4 string found at https://spatialreference.org/ref/epsg/nzgd2000-new-zealand-transverse-mercator-2000/
proj4str_nztm = '+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
schism_native = reader_netCDF_CF_unstructured_SCHISM.Reader(filename = '/data/schism_marl20080505_00z_3D.nc', proj4 = proj4str_nztm, use_3d = True)

# nordic_native = reader_ROMS_native.Reader(o.test_data_folder() +
#     '2Feb2016_Nordic_sigma_3d/Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc')
o.add_reader([reader_landmask,schism_native])

# Seed elements at defined positions, depth and time
o.seed_elements(lon=174.046669, lat=-40.928116, radius=0, number=10,
                z=np.linspace(0,-10, 10), time=schism_native.start_time)

#%%
# Running model
o.run(time_step=900, outfile='schism_native_output.nc')

#%%
# Print and plot results
print(o)
o.plot(linecolor='z', fast=True)
#o.animation()