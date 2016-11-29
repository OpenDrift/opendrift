#!/usr/bin/env python

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.pelagicegg import PelagicEggDrift

o = PelagicEggDrift(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=10, llcrnrlat=67.,
                    urcrnrlon=16, urcrnrlat=69.,
                    resolution='h', projection='merc')

o.add_reader([reader_norkyst, reader_basemap, reader_arome])

time = reader_arome.start_time

# spawn NEA cod eggs at defined position and time
o.seed_elements(14. , 68.1, z=-40, radius=2000, number=500,
                time=time, diameter=0.0014, neutral_buoyancy_salinity=31.25)
o.seed_elements(12.5, 68.5, z=-40, radius=2000, number=500,
                time=time, diameter=0.0014, neutral_buoyancy_salinity=31.25)


# Adjusting some configuration
o.config['processes']['turbulentmixing'] = True
o.config['turbulentmixing']['diffusivitymodel'] = 'windspeed_Sundby1983' # windspeed parameterization for eddy diffusivity
#o.config['turbulentmixing']['diffusivitymodel'] = 'environment' # use eddy diffusivity from ocean model 

# Vertical resolution and time step should be adjusted so to avoid getting
# output warnings like 'DEBUG: WARNING! some elements have p+q>1.'
o.config['turbulentmixing']['timestep'] = 10. # seconds
o.config['turbulentmixing']['verticalresolution'] = 2. # m
o.config['turbulentmixing']['TSprofiles'] = True # update T,S frequently in each inner loop for turbulent
o.config['turbulentmixing']['TSprofiles'] = False  # False is faster

# Running model (until end of driver data)
o.run(steps=6*2, time_step=1800)

# Print and plot results
print o

o.plot()
o.animation()
o.plot_vertical_distribution()

