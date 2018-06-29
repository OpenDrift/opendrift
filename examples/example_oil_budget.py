#!/usr/bin/env python

from datetime import timedelta

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil3D import OpenOil3D

o = OpenOil3D(loglevel=0)

# Arome
reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + 
#    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(
                    llcrnrlon=4, llcrnrlat=59.7,
                    urcrnrlon=7, urcrnrlat=61.5,
                    resolution='h', projection='merc')

#o.add_reader([reader_basemap, reader_norkyst, reader_arome])
o.fallback_values['x_wind'] = 7
o.fallback_values['y_wind'] = 0
o.fallback_values['x_sea_water_velocity'] = .7
o.fallback_values['y_sea_water_velocity'] = .3
#o.fallback_values['land_binary_mask'] = 0
#o.add_reader([reader_basemap, reader_norkyst])
o.add_reader([reader_basemap])

# Seeding some particles
lon = 4.8; lat = 60.0; # Outside Bergen

#time = [reader_arome.start_time,
#        reader_arome.start_time + timedelta(hours=30)]
time = reader_arome.start_time

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=3000, number=1000, time=time, z=0,
                oiltype='EKOFISK')

# Adjusting some configuration
o.set_config('processes:dispersion', True)
o.set_config('processes:evaporation', True)
o.set_config('processes:emulsification', True)
o.set_config('processes:turbulentmixing', True)
o.set_config('turbulentmixing:TSprofiles', False)
o.set_config('turbulentmixing:diffusivitymodel', 'windspeed_Sundby1983')
o.set_config('turbulentmixing:timestep', 2.) # seconds

# Running model (until end of driver data)
o.run(steps=4*20, time_step=900, export_buffer_length=10,
      outfile='oil_budget.nc')

# Print and plot results
print(o)
o.plot_oil_budget()
o.plot_property('water_fraction')
o.plot_property('water_fraction', mean=True)
o.plot()
o.plot_property('mass_oil')
o.plot_property('z')
o.plot_property('mass_evaporated')
o.plot_property('water_fraction')
o.plot_property('interfacial_area')
o.animation()
