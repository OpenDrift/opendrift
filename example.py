#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('/opdata_local/arome2_5/arome_metcoop_default2_5km_20150212_00.nc')

#reader_ec = reader_netCDF_CF_generic.Reader('/opdata/ec/la/nc4/ec_la_atmo_2015092200.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
#reader_norkyst = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')#, name='norkyst800_file')

# WAM10
#reader_wam10 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/wam/wam10_be')  # 6 hourly aggregates two months back
#reader_wam10 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/fou-hi/wam10h/g10kmwave.nc')

# Nordic4
#reader_nordic4 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/fou-hi/nordic4km-1h/Nordic-4km_SURF_1h_avg_00.nc')

# Arctic20
#reader_arctic20 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/arctic20km/1h/aggregate_be', name='arctic20_thredds')

# GlobCurrent
#reader_globcurrent = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREUL_HS-ALT_SUM-V01.0_FULL_TIME_SERIE')  # Total
#reader_globcurrent = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/GC_MOD_TIDE_GLO_010_FES2012_FULL_TIME_SERIE') # FES Tidal

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=2, llcrnrlat=59,
                    urcrnrlon=9, urcrnrlat=65,
                    resolution='h', projection='merc')
#reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=10.2, llcrnrlat=59.4,
#                    urcrnrlon=10.7, urcrnrlat=59.7,
#                    resolution='h', projection='merc')

#o.add_reader([reader_norkyst])
#o.add_reader([reader_arctic20, reader_arome, reader_basemap])
#o.add_reader([reader_norkyst, reader_arome, reader_basemap])
o.add_reader([reader_basemap, reader_norkyst, reader_arome])
#o.add_reader([reader_norkyst, reader_basemap])
#o.add_reader([reader_norkyst, reader_nordic4, reader_arctic20, reader_arome, reader_basemap])
#o.add_reader([reader_globcurrent, reader_basemap])

# Seeding some particles
#lon = 15; lat = 72.0; # Close to Norkyst boundary
#lon = 21; lat = 73.5; # Close to Norkyst boundary
#reader_norkyst.plot()
#lon = 10.6; lat = 57.33; # Laesoe, close to Norkyst boundary
#lon = 10.6; lat = 54.83; # outside Norkyst boundary
lon = 4.6; lat = 60.0; # Outside Bergen
#lon = 2.5; lat = 60.0; # Frigg/NOFO
lon = 6.73; lat = 62.78; # Outside Trondheim
#lat = 59.551010; lon = 10.618833 # Ytre Oslofjord
#lon = 10.546; lat = 59.486 # Godafoss
#lon = 9.76; lat = 58.94 # Full City
#lon = 9.76; lat = 58.84 # Full City south
#lon = 22.6; lat = 71.00; # Barents
#lon = 14.332435; lat = 68.045021; # Vestfjorden, Rohrs et al., 2012

time = None
#time = reader_wam10.start_time
#time = datetime(2015, 9, 22, 6, 0, 0)
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=4)]
#time = reader_arome.start_time

# Seed oil elements at defined position and time
o.seed_point(lon, lat, radius=0, number=1000, time=time)

print o

#o.seed_from_gml('/disk1/data/globoilrisk/ftp3.ksat.no/oilspills/RS2_20150608_171458_0045_SCNA_HH_SGF_401754_5438_11441747_Oil.gml', num_elements=1000)
#o.start_time = reader_nordic4.start_time

#o.set_projection(reader_arome.proj4)
#o.set_projection('+proj=latlong')

# Adjusting some configuration
o.config['drift']['wind_drift_factor'] = .01
o.config['processes']['diffusion'] = True
o.config['processes']['dispersion'] = True
o.config['processes']['evaporation'] = True
o.config['processes']['emulsification'] = True
o.config['drift']['current_uncertainty'] = .0
o.config['drift']['wind_uncertainty'] = 2

# Running model (until end of driver data)
#o.run(steps=80*4, time_step=900)
o.run(steps=80*4, time_step=900, outfile='openoil.nc')

print o.history

# Print and plot results
print o
#o.plot(background=['x_sea_water_velocity', 'y_sea_water_velocity'], buffer=.5)
o.plot()
#o.plot_property('mass_oil')
#o.plot_property('x_sea_water_velocity')
