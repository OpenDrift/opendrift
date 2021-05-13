#!/usr/bin/env python
"""
Droplet distribution
==================================

Plotting different droplet size distributions used in Opendrift
"""

import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime, timedelta
from opendrift.models.openoil import OpenOil

#%%
# droplet size interval for plotting
dmin = 1.e-6
dmax = 3.e-3

#%%
# simulation with Johansen et al. (2015) droplet spectrum
o = OpenOil(loglevel=20, weathering_model='noaa')
o.set_config('environment:fallback:land_binary_mask', 0)
o.set_config('environment:fallback:x_sea_water_velocity', -.2)
o.set_config('environment:fallback:y_sea_water_velocity', 0)
o.set_config('environment:fallback:x_wind', 6.)
o.set_config('environment:fallback:y_wind', 0)
o.set_config('environment:fallback:sea_water_temperature', 5.)
o.set_config('environment:fallback:sea_surface_wave_stokes_drift_x_velocity', .3)
o.set_config('environment:fallback:sea_surface_wave_stokes_drift_y_velocity', 0)
o.set_config('wave_entrainment:droplet_size_distribution', 'Johansen et al. (2015)')
o.set_config('processes:evaporation', False)
o.set_config('processes:dispersion', False)
o.seed_elements(lon=4, lat=60, time=datetime.utcnow(), number=10000, radius=100,
                 z=0, oil_type='TROLL, STATOIL', oil_film_thickness=0.005)
o.run(duration=timedelta(hours=2), time_step=3600)

droplet_diameters = o.elements.diameter
sd = 0.4
Sd = np.log(10.)*sd
DV_50 = np.median(droplet_diameters)
DN_50 = np.exp( np.log(DV_50) - 3*Sd**2 )

print('DV_50: %f' % DV_50)
print('DN_50: %f' % DN_50)

#%%
# Plotting
plt.figure(figsize=[14,14])
plt.subplot(3, 2, 1)
nVpdf, binsV, patches = plt.hist(droplet_diameters, 100, range=(dmin,dmax), align='mid')
plt.xlabel('Droplet diameter d [m]', fontsize=8)
plt.ylabel('V(d)', fontsize=8)
plt.title('volume spectrum\nJohansen et al. (2015) distribution in OpenOil', fontsize=10)

plt.subplot(3,2,2)
nVcum, binsV, patches = plt.hist(droplet_diameters, 100, range=(dmin,dmax), align='mid', cumulative=True)
plt.xlabel('Droplet diameter d [m]', fontsize=8)
plt.ylabel('V(d)', fontsize=8)
plt.title('cumulative volume spectrum', fontsize=10)

#%%
# calculate number spectrum from volume spectrum
d = 0.5* (binsV[1:] + binsV[:-1])
V = 4./3. * np.pi * (d/2.)**3
Npdf = nVpdf / V
Ncum = np.cumsum(Npdf)

#%%
# calculate theoretical cumulative Number distribution
spectrum = (np.exp(-(np.log(d) - np.log(DN_50))**2 / (2 * Sd**2))) / (d * Sd * np.sqrt(2 * np.pi)) # from OpenOil median diameter
DN_50_johansen = 0.000483
spectrumJ = (np.exp(-(np.log(d) - np.log(DN_50_johansen))**2 / (2 * Sd**2))) / (d * Sd * np.sqrt(2 * np.pi)) # from parameters in Johansen et al. 2015

spectrum_cum = np.cumsum(spectrum)
spectrumJ_cum = np.cumsum(spectrumJ)
# calculate theoretical number distribution pdf
#spectrum_pdf = spectrum / np.sum(spectrum)
#spectrumJ_pdf = spectrumJ / np.sum(spectrumJ)


plt.subplot(3,2,3)
plt.plot(d, Npdf/np.sum(Npdf), lw=2)
#plt.plot(d, spectrum/np.sum(spectrum), label='curve fit from OpenOil', lw=2)
plt.plot(d, spectrumJ/np.sum(spectrumJ), label='curve fit from Johansen et al. 2015', lw=2)
plt.xlabel('Droplet diameter d [m]', fontsize=8)
plt.ylabel('N(d)', fontsize=8)
plt.title('number spectrum', fontsize=10)

plt.subplot(3,2,4)
plt.plot(d, Ncum/np.max(Ncum), label='OpenOil result', lw=2)
#plt.plot(d, spectrum_cum/np.max(spectrum_cum), label='OpenOil par.', lw=2)
plt.plot(d, spectrumJ_cum/np.max(spectrumJ_cum), label='Johansen et al. 2015', lw=2)
plt.xlabel('Droplet diameter d [m]', fontsize=8)
plt.ylabel('N(d)', fontsize=8)
plt.title('cumulative number spectrum', fontsize=10)
plt.legend(loc='lower right')

plt.subplot(3,2,5)
plt.loglog(d, Npdf/np.sum(Npdf), lw=2)
#plt.loglog(d, spectrum/np.sum(spectrum), label='curve fit from OpenOil', lw=2)
plt.loglog(d, spectrumJ/np.sum(spectrumJ), label='curve fit from Johansen et al. 2015', lw=2)
plt.xlabel('Droplet diameter d [m]', fontsize=8)
plt.ylabel('N(d)', fontsize=8)
plt.title('number spectrum', fontsize=10)

plt.subplot(3,2,6)
plt.semilogx(d, Ncum/np.max(Ncum), label='OpenOil result', lw=2)
#plt.semilogx(d, spectrum_cum/np.max(spectrum_cum), label='OpenOil par.', lw=2)
plt.semilogx(d, spectrumJ_cum/np.max(spectrumJ_cum), label='Johansen et al. 2015', lw=2)
plt.xlabel('Droplet diameter d [m]', fontsize=8)
plt.ylabel('N(d)', fontsize=8)
plt.title('cumulative number spectrum', fontsize=10)
#plt.legend(loc='lower right')
plt.show()



