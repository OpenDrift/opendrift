#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil3D import OpenOil3D

################################################################
# NOTE: This example requires some huge 3D input files, 
#       which are not included in OpenDrift repository. 
#       Contact Knut-Frode (knutfd@met.no) if you would like
#       access the data to run this example (~5 GB)
################################################################

o = OpenOil3D(loglevel=0)  # Set loglevel to 0 for debug information

ncfile = 'oil3Dmixing.nc'
import_file = False

if import_file is True:
    o.io_import_file(ncfile)
else:
    # Arome
    reader_arome = reader_netCDF_CF_generic.Reader('/disk2/data/opendrift_test_data/15jan2016/AROME_MetCoOp_00_DEF.nc')
    #reader_arome = reader_netCDF_CF_generic.Reader('/opdata/arome2_5_main/AROME_MetCoOp_00_DEF.nc')

    # Norkyst
    reader_norkyst = reader_netCDF_CF_generic.Reader('/disk2/data/opendrift_test_data/15jan2016/NorKyst-800m_ZDEPTHS_his_00.nc')
    #reader_norkyst = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')

    # Landmask (Basemap)
    # Stad
    reader_basemap = reader_basemap_landmask.Reader(
                        llcrnrlon=3.0, llcrnrlat=61.4,
                        urcrnrlon=6.2, urcrnrlat=62.6,
                        resolution='h', projection='merc')

    o.add_reader([reader_norkyst, reader_basemap, reader_arome])

    # Experiment with constant wind instead of arome
    #o.add_reader([reader_norkyst, reader_basemap])
    #o.fallback_values['x_wind'] = 20
    #o.fallback_values['y_wind'] = 20

    # Seeding some particles
    #lon = 5.0; lat = 60.0; # Outside Bergen
    #lon = 5.15; lat = 58.85; # Outside Stavanger
    lon = 4.9; lat = 62.1; # Stad
    #lon = 13.8; lat = 67.57; # Vestfjorden

    #time = [reader_arome.start_time,
    #        reader_arome.start_time + timedelta(hours=30)]
    time = reader_arome.start_time

    # Seed oil elements at defined position and time
    o.seed_elements(lon, lat, z=-0.5, radius=0, number=3000, time=time,
                    density=880)

    # Adjusting some configuration
    o.config['drift']['wind_drift_factor'] = .02

    o.config['processes']['turbulentmixing'] = True
    o.config['turbulentmixing']['diffusivitymodel'] = 'windspeed_Sundby1983'
    #o.config['turbulentmixing']['diffusivitymodel'] = 'stepfunction'
    o.config['turbulentmixing']['timestep'] = 2. # seconds

    o.config['processes']['diffusion'] = False
    o.config['drift']['current_uncertainty'] = .1
    o.config['drift']['wind_uncertainty'] = 2
    o.config['processes']['evaporation'] = True
    o.config['processes']['dispersion'] = True
    o.config['processes']['emulsification'] = True

    # Running model (until end of driver data)
    o.run(steps=66*2, time_step=1800, outfile=ncfile)

###########################
# Print and plot results
###########################
print o

o.plot_oil_budget()
o.plot()
#o.plot(linecolor='z')
o.animation()
#o.plot_vertical_distribution()
o.plot_vertical_distribution_slider()

