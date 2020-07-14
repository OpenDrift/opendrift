"""
3D - Derivadores
=============
"""
#!/usr/bin/env python
import numpy as np

from opendrift.readers import reader_ECOM_3D
from opendrift.models.oceandrift3D import OceanDrift3D

o = OceanDrift3D(loglevel=0)	  # Set loglevel to 0 for debug information

reader_pcse = reader_ECOM_3D.Reader('/home/arian/IC/Deriva/Arian_cases/IC_drifters/ECOM_data/Output_netcdf/run262_in_z.nc')

o.add_reader([reader_pcse])

#DEVS COMEÇANDO NO PRIMEIRO PONTO DE PART_LOCATION

# Dev 1
#lon = -45.21035; lat = -26.01430;
#z= -1.43061


#Dev 2
#lon = -45.14352 ; lat = -26.08981 ;
#z =  -1.69214

#Dev 3
#lon = -45.06477 ; lat =-26.1738275  ;
#z = -1.91010

#Dev 4
#lon =-44.97670   ; lat = -26.2615425 ;
#z = -2.09853
	
#Dev 5
#lon = -44.88531  ; lat = -26.34655  ;
#z = -2.55560

#Dev 6
#lon = -45.58546  ; lat = -26.21616 ;
#z = -1.86552

#Dev 7
#lon = -45.29894  ; lat = -26.66503 ;
#z = -2.56940

#Dev 8
#lon =  -46.54288 ; lat = -27.02260 ;
#z = -2.59892	

#Dev 9
lon =  -46.01746  ; lat = -27.85090  ;
z =-2.56099

time = reader_pcse.start_time
# Seed elements at defined position and time


o.seed_elements(lon, lat, z = z, radius=1, number=2, time=time)

# Adjusting some configuration
#o.set_config('drift:current_uncertainty', .1)
#o.set_config('drift:wind_uncertainty', 2)
#o.set_config('processes:verticaladvection', True)
#o.set_config('process:turbulentmixing', True)

# Running model (until end of driver data)
#o.run(duration=timedelta(hours=24))
o.run(steps=282, time_step=900, time_step_output=900, outfile='/home/arian/IC/arquivos_cdf/Dev262_9_noangle_nowind_Opendrift.nc') 
#o.run(steps=283, time_step=900)   
#o.run(duration=timedelta(hours=24))
#o.run(steps = 235, time_step = 1800)

# Print and plot results
#print(o)
o.plot(linecolor='z', cmap = 'brg', filename = 'Dev262_9')  # Color lines according to depth
#o.plot_property('z')

#Jet colors with angle
print(o)
#o.plot(linecolor='z', cmap = 'jet')  # Color lines according to depth
#o.plot_property('z')




#
