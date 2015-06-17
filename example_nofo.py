#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

o = OpenOil(loglevel=0)  # Set loglevel to 0 for debug information

# Arome
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
reader_arome1 = reader_netCDF_CF_generic.Reader('/vol/hindcast3/nofo_olje_paa_vann_2015/arome_metcoop_default2_5km_20150610_00.nc')
reader_arome2 = reader_netCDF_CF_generic.Reader('/vol/hindcast3/nofo_olje_paa_vann_2015/arome_metcoop_default2_5km_20150610_06.nc')
reader_arome3 = reader_netCDF_CF_generic.Reader('/vol/hindcast3/nofo_olje_paa_vann_2015/arome_metcoop_default2_5km_20150611_06.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('/opdata_local/arome2_5/arome_metcoop_default2_5km_20150212_00.nc')

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
#reader_norkyst = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')#, name='norkyst800_file')

# WAM10
#reader_wam10 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/wam/wam10_be')  # 6 hourly aggregates two months back
reader_wam10 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/fou-hi/wam10h/g10kmwave.nc')

# Nordic4
reader_nordic4 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/fou-hi/nordic4km-1h/Nordic-4km_SURF_1h_avg_00.nc')

# Arctic20
#reader_arctic20 = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/arctic20km/1h/aggregate_be', name='arctic20_thredds')

# GlobCurrent
#reader_globcurrent = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/CLS-L4-CUREUL_HS-ALT_SUM-V01.0_FULL_TIME_SERIE')  # Total
#reader_globcurrent = reader_netCDF_CF_generic.Reader('http://tds0.ifremer.fr/thredds/dodsC/GC_MOD_TIDE_GLO_010_FES2012_FULL_TIME_SERIE') # FES Tidal

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=-5, llcrnrlat=54,
                    urcrnrlon=27, urcrnrlat=70, resolution='h')
#reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=0, llcrnrlat=54,
#                    urcrnrlon=6, urcrnrlat=59, resolution='h')

#o.add_reader([reader_norkyst])
#o.add_reader([reader_arctic20, reader_arome, reader_basemap])
o.add_reader([reader_norkyst, reader_arome3, reader_arome2, reader_arome1, reader_basemap])
#o.add_reader([reader_norkyst, reader_basemap])
#o.add_reader([reader_arctic20, reader_basemap])
#o.add_reader([reader_norkyst, reader_nordic4, reader_arctic20, reader_arome, reader_basemap])
#o.add_reader([reader_norkyst, reader_arome, reader_basemap])

# Seeding some particles
#lon = 10.546; lat = 59.486 # Godafoss
#lon = 9.76; lat = 58.94 # Full City
#lon = 9.76; lat = 58.84 # Full City south
#lon = 22.6; lat = 71.00; # Barents
#lon = 14.332435; lat = 68.045021; # Vestfjorden, Rohrs et al., 2012
lon = 4.5; lat = 60  # Bergen
time = None
#time = reader_wam10.start_time
#time = datetime(2015, 6, 9, 9, 0, 0)
#time = datetime(2011, 7, 1, 0, 0, 0)
#o.seed_point(lon, lat, radius=100000, number=1000, time=time)
import numpy as np
kwargs = {}
#kwargs['lon'] = np.array([2.385876, 2.391800])
o.seed_point(2.385876, 60.032848, radius=0, number=500, time=time)
#kwargs['lat'] = np.array([60.032848, 60.051200])
o.seed_point(2.391800, 60.051200, radius=0, number=500, time=time)
#kwargs['ID'] = np.array([1, 2])
#kwargs['age_exposure_seconds'] = np.array([6, 6])*3600  # Some already avaporated
#kwargs['fraction_evaporated'] = np.array([0.3, 0.3])  # Some already avaporated
#o.elements = o.ElementType(**kwargs)
o.start_time = datetime(2015, 6, 10, 5, 15, 0)

#o.seed_from_gml('/disk1/data/globoilrisk/ftp3.ksat.no/oilspills/RS2_20150611_172613_0086_SCNB_HH_SGF_402348_0016_11441894_Oil.gml', num_elements=1000)
o.seed_from_gml('/disk1/data/globoilrisk/ftp3.ksat.no/oilspills/RS2_20150615_053843_0045_SCNA_HHHV_SGF_403055_5578_11441860_Oil.gml', num_elements=1000)
#o.start_time = reader_nordic4.start_time

#o.set_projection(reader_arome.proj4)
#o.set_projection('+proj=latlong')

# Running model (until end of driver data)
o.wind_factor = 0.035
o.run(steps=4*10, time_step=900, outfile='openoil.nc')

# Print and plot results
print o
#o.plot(background='x_sea_water_velocity', buffer=.5)
#o.plot(background=['x_sea_water_velocity', 'y_sea_water_velocity'], buffer=.5)
#dfiles = ['/disk1/data/nofo_olje_pa_vann_2015/ovp_isphere1-300234062952180-20150612T090517UTC.csv', '/disk1/data/nofo_olje_pa_vann_2015/ovp_isphere2-300234062957180-20150612T090518UTC.csv']
#dfiles = ['/disk1/data/nofo_olje_pa_vann_2015/ovp_isldmb1-300234062953170-20150612T090518UTC.csv', '/disk1/data/nofo_olje_pa_vann_2015/ovp_isldmb2-300234062959190-20150612T090518UTC.csv']
#dfiles = ['/disk1/data/nofo_olje_pa_vann_2015/isphere1.csv', '/disk1/data/nofo_olje_pa_vann_2015/isphere2.csv']
#o.plot(buffer=.05, drifter_file=dfiles)
o.plot(buffer=.05)
o.plot_property('fraction_evaporated')
o.plot_property('water_content')
