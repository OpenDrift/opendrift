#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
reader_arome = reader_netCDF_CF_generic.Reader('test_data/16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst = reader_netCDF_CF_generic.Reader('test_data/16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

# WAM10
#reader_wam10 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/wam/wam10_be')  # 6 hourly aggregates two months back
#reader_wam10 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/fou-hi/wam10h/g10kmwave.nc')

# GlobCurrent
#reader_globcurrent = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREUL_HS-ALT_SUM-V01.0_FULL_TIME_SERIE')  # Total
#reader_globcurrent = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/GC_MOD_TIDE_GLO_010_FES2012_FULL_TIME_SERIE') # FES Tidal

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=4, llcrnrlat=59.8,
                    urcrnrlon=6, urcrnrlat=61,
                    resolution='h', projection='merc')

o.add_reader([reader_basemap, reader_norkyst, reader_arome])

# Seeding some particles
lon = 4.6; lat = 60.0; # Outside Bergen
#lon = 6.73; lat = 62.78; # Outside Trondheim

time = None
#time = datetime(2015, 9, 22, 6, 0, 0)
time = [reader_arome.start_time,
        reader_arome.start_time + timedelta(hours=30)]
#time = reader_arome.start_time

# Seed oil elements at defined position and time
o.seed_elements(lon, lat, radius=50, number=3000, time=time)

print o

#o.seed_from_gml('/disk1/data/globoilrisk/ftp3.ksat.no/oilspills/RS2_20150608_171458_0045_SCNA_HH_SGF_401754_5438_11441747_Oil.gml', num_elements=1000)
#o.start_time = reader_nordic4.start_time

#o.set_projection(reader_arome.proj4)
#o.set_projection('+proj=latlong')

# Adjusting some configuration
o.config['drift']['wind_drift_factor'] = .02
o.config['processes']['diffusion'] = True
o.config['processes']['dispersion'] = True
o.config['processes']['evaporation'] = True
o.config['processes']['emulsification'] = True
o.config['drift']['current_uncertainty'] = .1
o.config['drift']['wind_uncertainty'] = 2

# Running model (until end of driver data)
o.run(steps=66*2, time_step=1800, outfile='openoil.nc')

# Print and plot results
print o
#o.plot(background=['x_sea_water_velocity', 'y_sea_water_velocity'], buffer=.5)
o.animation()
#o.animation(filename='openoil_time_seed.gif')
o.plot()
#o.plot_property('mass_oil')
#o.plot_property('x_sea_water_velocity')
