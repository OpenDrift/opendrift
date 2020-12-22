#!/usr/bin/env python
"""
LCS Norkyst
==================================
"""

from datetime import datetime, timedelta
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information

# Norkyst ocean model
reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')

#o.add_reader([reader_norkyst, reader_arome])
o.add_reader([reader_norkyst])


#%%

x0,y0 = reader_norkyst.lonlat2xy(4.5,60.5)
d = 30000
lcsproj = reader_norkyst.proj
lcs = o.calculate_lcs(
    time=reader_norkyst.start_time + timedelta(hours=24),
    time_step=timedelta(minutes=15), reader=lcsproj,
    duration=timedelta(hours=3), delta=400, domain=[x0-d,x0+d,y0-d,y0+d])

l1 = lcs['eigval'][0,:,:,0]
l2 = lcs['eigval'][0,:,:,1]
mask = l2>l1

lmax = l1.copy()
lmax[mask] = l2[mask]

lmin = l2.copy()
lmin[mask] = l1[mask]

xi1 = lcs['eigvec'][0,:,:,0]
xi2 = lcs['eigvec'][0,:,:,1]

stretch = xi1.copy()
stretch[mask] = xi2[mask]

shrink = xi2.copy()
shrink[mask] = xi1[mask]


# find local maxima of largest eigenvalue
from skimage.feature import peak_local_max
#peaks = peak_local_max(lcs['eigval'][0,:,:,1],5)
peaks = peak_local_max(-lmin,5)
plon = [lcs['lon'][peaks[i,0],peaks[i,1]] for i in range(peaks.shape[0])]
plat = [lcs['lat'][peaks[i,0],peaks[i,1]] for i in range(peaks.shape[0])]


from opendrift.readers.reader_constant_2d import Reader
x,y = reader_norkyst.lonlat2xy(lcs['lon'],lcs['lat'])
r = Reader(x=x, y=y, proj4=reader_norkyst.proj4, array_dict = {'x_sea_water_velocity': stretch[:,:,0], 'y_sea_water_velocity': stretch[:,:,1]})


xi = OceanDrift(loglevel=20)
xi.add_reader(r)
xi.seed_elements(plon, plat, time=datetime.now())
xi.run(duration=timedelta(hours=3), time_step=300)
str_lon = xi.history['lon']
str_lat = xi.history['lat']

#xi.plot(linewidth=2, linecolor='r',show_particles=False)

xi.reset()
xi.seed_elements(plon, plat, time=datetime.now())
xi.run(duration=timedelta(hours=3), time_step=-300)
str_lon = np.concatenate((str_lon, xi.history['lon'][:,::-1]), axis=0)
str_lat = np.concatenate((str_lat, xi.history['lat'][:,::-1]), axis=0)


#%%
# Simulation from beginning and up to 30 hours (time of LCS)

o.reset()
lons = np.linspace(4.0, 5.0, 100)
lats = np.linspace(60., 61., 100)
lons, lats = np.meshgrid(lons, lats)
lons = lons.ravel()
lats = lats.ravel()
o.seed_elements(lons, lats, radius=0, number=10000,
                time=reader_norkyst.start_time)

o.run(end_time=reader_norkyst.start_time+timedelta(hours=24),
      time_step=timedelta(minutes=30))
ax, plt = o.plot(cmap='Reds', vmax=2, markersize=1, colorbar=True, show_particles=True, show_initial=False, linewidth=0, show=False)

'''
fig = plt.figure()
#proj=reader_norkyst.proj

lon, lat = lcs['lon'], lcs['lat']
x,y = reader_norkyst.lonlat2xy(lon, lat)
px,py = reader_norkyst.lonlat2xy(plon, plat)
#strx, stry = reader_norkyst.lonlat2xy(str_lon, str_lat)
fig.add_subplot(121)
plt.pcolormesh(x, y, np.log(np.sqrt(lmin)),vmin=-1,vmax=3, cmap='cividis'),plt.colorbar()
plt.quiver(x[::5,::5], y[::5,::5], shrink[::5,::5,0], shrink[::5,::5,1])
plt.plot(px, py, '.r')
plt.title('shrink lines')
fig.add_subplot(122)

plt.pcolormesh(x,y , np.log(np.sqrt(lmax)), cmap='cividis'),plt.colorbar()
plt.quiver(x[::5,::5], y[::5,::5], stretch[::5,::5,0], stretch[::5,::5,1])

plt.title('stretch lines')

plt.figure()
'''
gcrs = ccrs.PlateCarree()
ax.plot(str_lon.T, str_lat.T, '-r', ms=0.5, transform=gcrs)
plt.show()

#o.animation(buffer=0, lcs=ftle, hide_landmask=True)
