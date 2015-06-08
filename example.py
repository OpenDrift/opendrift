#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('/opdata_local/arome2_5/arome_metcoop_default2_5km_20150212_00.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
#reader_norkyst = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')#, name='norkyst800_file')

# WAM10
#reader_wam10 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/wam/wam10_be')  # 6 hourly aggregates two months back
#reader_wam10 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/fou-hi/wam10h/g10kmwave.nc')

# Nordic4
reader_nordic4 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/fou-hi/nordic4km-1h/Nordic-4km_SURF_1h_avg_00.nc')

# Arctic20
reader_arctic20 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/arctic20km/1h/aggregate_be', name='arctic20_thredds')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=-5, llcrnrlat=54,
                    urcrnrlon=27, urcrnrlat=79, resolution='h')

#o.add_reader([reader_norkyst])
#o.add_reader([reader_arctic20, reader_arome, reader_basemap])
#o.add_reader([reader_norkyst, reader_arome, reader_basemap])
#o.add_reader([reader_norkyst, reader_basemap])
#o.add_reader([reader_arctic20, reader_basemap])
#o.add_reader([reader_norkyst, reader_arctic20, reader_arome, reader_wam10, reader_basemap])
o.add_reader([reader_norkyst, reader_nordic4, reader_arctic20, reader_arome, reader_basemap])

# Seeding some particles
#lon = 15; lat = 72.0; # Close to Norkyst boundary
#lon = 21; lat = 73.5; # Close to Norkyst boundary
#reader_norkyst.plot()
#lon = 10.6; lat = 57.33; # Laesoe, close to Norkyst boundary
#lon = 10.6; lat = 54.83; # outside Norkyst boundary
lon = 4.9; lat = 59.9; # Outside Bergen
#lon = 2.0; lat = 60.0; # Frigg/NOFO
#lon = 6.73; lat = 62.78; # Outside Trondheim
#lon = 10.546; lat = 59.486 # Godafoss
#lon = 9.76; lat = 58.94 # Full City
#lon = 9.76; lat = 58.84 # Full City south
#lon = 22.6; lat = 71.00; # Barents
time = None
#time = reader_wam10.start_time
#time = datetime(2015, 6, 9, 12, 0, 0)
o.seed_point(lon, lat, radius=10000, number=1000, massOil=5, time=time)

#o.seed_from_gml('/disk1/data/globoilrisk/ftp3.ksat.no/oilspills/RS2_20150604_173047_0076_SCWA_VV_SGF_400947_5309_11313456_Oil.gml', num_elements=1000)
#o.seed_from_gml('/disk1/data/globoilrisk/ftp3.ksat.no/oilspills/RS2_20150531_061945_0086_SCNA_HHHV_SGF_400108_5164_11293542_Oil.gml', num_elements=1000)
#o.start_time = reader_nordic4.start_time

#o.set_projection(reader_arome.proj4)
#o.set_projection('+proj=latlong')

# Running model (until end of driver data)
o.run(steps=80, time_step=1800, outfile='openoil.nc')

# Print and plot results
print o
#o.plot(background='sea_water_potential_temperature')
o.plot(buffer=.5)
