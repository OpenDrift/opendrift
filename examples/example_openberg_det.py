#!/usr/bin/env python

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_basemap_landmask

from opendrift.models.openberg import OpenBerg



######################
# Initialize model
######################
steps = 60  # This is the number of forecast steps

o = OpenBerg()  # Basic drift model suitable for icebergs

############################################
# Preparing Readers 
############################################

reader_current = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
   '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

reader_wind = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
   '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
   
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=2.9, llcrnrlat=59.7,
                        urcrnrlon=4.9, urcrnrlat=61.5, resolution='l',
                        projection='gall')    
o.add_reader([reader_current,reader_wind,reader_basemap])

#######################
# Seeding elements
#######################

# Icebergs are moved with the ocean current as per Barker et al (2004), 
# in addition to a fraction of the wind speed (wind_drift_factor). 
# This factor depends on the properties of the elements. 
# Default empirical values are:
# - Wind drift fraction: 0.018 (1.8 %) (Garret 1985)
# - Iceberg size: 	Keel dept = 60m
#					Waterline length = 90.5m


o.seed_elements(3.3, 60., radius=3000,
                time=reader_current.start_time,water_line_length=100,
                keel_depth=90,number=100)

                   
print('Starting free run .../n')
    
print('Start time: ' + str(o.start_time)) 

#######################
# Running model
#######################    

o.run(time_step=3600,steps=steps)

#########################    
# Print and plot results
#########################
o.plot(filename='example_det.pdf')
o.animation(filename='example_det.mp4')

print('############## Latitudes:', o.history['lat'])
print('############## Longitudes:',o.history['lon']) 

 
    
  



