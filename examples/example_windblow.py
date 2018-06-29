#!/usr/bin/env python

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.windblow import WindBlow

o = WindBlow(loglevel=0)  # Set loglevel to 0 for debug information

##############################################################
# Example of elements blowing with the wind, also over land
##############################################################

#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc')
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '2Feb2016_Nordic_sigma_3d/AROME_MetCoOp_00_DEF.nc_20160202_subset')

o.add_reader([reader_arome])

# Seeding some particles
lat = 68.5; lon = 16.0  # Lofoten
o.seed_elements(lon, lat, radius=5000, number=1000,
                time=reader_arome.start_time)

# Running model
o.run(steps=48*4, time_step=900)

# Print and plot results
print(o)
o.plot(buffer=.5)
o.animation()
