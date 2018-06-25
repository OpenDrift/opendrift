#!/usr/bin/env python

from datetime import datetime

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

o.add_reader([reader_norkyst, reader_arome])
o.set_config('processes:evaporation', False)

############################################################
# Seed oil particles within contour detected from satellite
############################################################
o.seed_from_gml(o.test_data_folder() + 'radarsat_oil_satellite_observation/RS2_20151116_002619_0127_SCNB_HH_SGF_433012_9730_12182143_Oil.gml', num_elements=2000)
############################################################
# Additional continous point release, lasting 24 hours
############################################################
o.seed_elements(3.8, 60.9, radius=0, number=1000,
                time=[datetime(2015,11,16,8), datetime(2015,11,17,8)]) 
############################################################
# Additional cone release (e.g. from moving ship)
############################################################
o.seed_elements([3.6, 4.4], [61.5, 61.2], radius=[1000, 10000],
                number=1000, cone=True,
                time=[datetime(2015,11,16,1), datetime(2015,11,16,8)])

# Running model (until end of driver data)
o.run(steps=50*4, time_step=900)

# Print and plot results
print(o)
o.plot()
o.animation()
