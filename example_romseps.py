#!/usr/bin/env python

from datetime import datetime, timedelta

from readers import reader_basemap_landmask
from readers import reader_netCDF_CF_generic
from models.openoil import OpenOil

import numpy as np
import multiprocessing
#nproc=multiprocessing.cpu_count() #number of parallel threads
nproc=1 #number of parallel threads
print 'Will use '+str(nproc)+' processors if run in parallel'
pool = multiprocessing.Pool(processes=nproc)

def write_ascii(object, outfile, outpath='./'):
    # Write to file:
    for t in range(len(object.lats[:,0])):
        f = open(outpath+'tstep_'+str(t)+'_'+outfile, 'w')
        f.write('%% particleID,lat,lon\n')
        for p in range(len(object.lats[0,:])):
            f.write(str(p)+', %3.2f, %3.2f\n' % (object.lats[t,p],object.lons[t,p]))
        f.close()


numeps = 11
#time = datetime(2015, 3, 17, 12, 0)
time = datetime.now()
lon = 10.0; lat = 58.0; # Skagerrak
#o = np.empty(numeps,dtype=object)


# Arome
reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('/opdata/arome2_5_main/AROME_MetCoOp_00_DEF.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
#reader_norkyst = reader_netCDF_CF_generic.Reader('/opdata/roms/NorKyst-800m_ZDEPTHS_his_00.nc')#, name='norkyst800_file')

# Landmask (Basemap)
reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=-5, llcrnrlat=54,
                    urcrnrlon=20, urcrnrlat=79, resolution='i')

def do_eps_drift(i):
    o = OpenOil(loglevel=20)  # Set loglevel to 0 for debug information
    reader_eps = reader_netCDF_CF_generic.Reader('/opdata/romseps/ocean_avg_z_'+str(i)+'.nc')
    o.add_reader([reader_arome, reader_eps, reader_basemap])
    #print o
    #time = None
    o.seed_elements(lon, lat, radius=10000, number=100, massOil=5, time=time)
    o.run(steps=120)
    #print o
    print 'done number '+str(i)
    return o


#o = pool.map(do_eps_drift,range(numeps))
o = map(do_eps_drift,range(numeps))


#print o

for i in range(numeps):
    #o[i].plot('eps_member'+str(i)+'.png',latmin=57.5,latmax=58.7,lonmin=9.5,lonmax=11)
    write_ascii(o[i],'member'+str(i)+'.dat', outpath='eps_out/')
    print 'Printed #'+str(i)

