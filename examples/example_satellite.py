#!/usr/bin/env python

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=3.5, llcrnrlat=60.4,
                    urcrnrlon=7, urcrnrlat=61.8,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_norkyst, reader_arome])

############################################################
# Seed oil particles within contour detected from satellite
############################################################
o.seed_from_gml(o.test_data_folder() + 'radarsat_oil_satellite_observation/RS2_20151116_002619_0127_SCNB_HH_SGF_433012_9730_12182143_Oil.gml',
    num_elements=2000)

# Adjusting some configuration
o.set_config('processes:dispersion', True)
o.set_config('processes:evaporation', False)
o.set_config('processes:emulsification', True)
o.set_config('drift:current_uncertainty', .3)  # Diffusion
o.set_config('drift:wind_uncertainty', 2)

# Running model (until end of driver data)
o.run(steps=66*4, time_step=900, outfile='openoil.nc')

# Print and plot results
print(o)
o.animation()
o.plot()
