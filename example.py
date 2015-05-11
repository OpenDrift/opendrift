#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('/opdata_local/arome2_5/arome_metcoop_default2_5km_20150212_00.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
#reader_norkyst = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')#, name='norkyst800_file')

# Arctic20
reader_arctic20 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/arctic20km/1h/aggregate_be', name='arctic20_thredds')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=-5, llcrnrlat=54,
                    urcrnrlon=27, urcrnrlat=79, resolution='h')

#o.add_reader([reader_norkyst])
#o.add_reader([reader_arctic20, reader_arome, reader_basemap])
#o.add_reader([reader_norkyst, reader_basemap])
#o.add_reader([reader_arctic20, reader_basemap])
#o.add_reader([reader_norkyst, reader_arome, reader_basemap])
o.add_reader([reader_norkyst, reader_arctic20, reader_arome, reader_basemap])

# Seeding some particles
#lon = 15; lat = 72.0; # Close to Norkyst boundary
lon = 21; lat = 73.5; # Close to Norkyst boundary
#reader_norkyst.plot()
#lon = 10.6; lat = 57.33; # Laesoe, close to Norkyst boundary
#lon = 10.6; lat = 54.83; # outside Norkyst boundary
lon = 4.9; lat = 60.0; # Outside Bergen
#lon = 6.73; lat = 62.78; # Outside Trondheim
#lon = 22.6; lat = 71.00; # Barents
time = None
#time = datetime(2015, 5, 11, 0, 0, 0)
#time = reader_norkyst.start_time
o.seed_point(lon, lat, radius=5000, number=5, massOil=5, time=time)

print o
#o.set_projection('+proj=longlat')

# Running model (until end of driver data)
o.run(steps=250, time_step=900)

# Print and plot results
print o
#o.plot(background='sea_water_potential_temperature')
o.plot(buffer=.5)
