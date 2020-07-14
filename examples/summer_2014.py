"""
3D - Derivadores
=============
"""
#!/usr/bin/env python
import numpy as np

from opendrift.readers import reader_ECOM_3D
from opendrift.models.oceandrift3D import OceanDrift3D

o = OceanDrift3D(loglevel=0)  # Set loglevel to 0 for debug information

reader_pcse = reader_ECOM_3D.Reader('/home/arian/IC/arquivos_cdf/summer2014_in_z.nc')
#reader_pcse = reader_ECOM_3D.Reader('/home/arian/IC/arquivos_cdf/run261_in_z.nc')
o.add_reader([reader_pcse])

#Posições


#Santos
#lat = -24.05   ; lon = -46.4 ;
#z = -15

#São Sebastião
lat_ss = -23.91 ; lon_ss = -45.71;
z_ss = -15

#Ilha Grande
lat_ig = - 23.56 ; lon_ig = -45.035;
z_ig = -15

#Cananeia
lat_cn = -25.06 ; lon_cn = -47.832;
z_cn = -15

#Offshore

#Off1
#lat_1 = -25.24 ; lon_1 = -47.1;
#z_1 = -15

#Off2
lat_2 = -25.22 ; lon_2 = -44.6;
z_2 = -15

#Off3
#lat_3 = -24.73 ; lon_3 = -43;
#z_3 = -15


#Off5
lat_5 = -23.05 ; lon_5 = -40.25;
z_5 = -15




time_start = reader_pcse.start_time
time_delay = reader_pcse.delay_time
# Seed elements at defined position and time


#o.seed_elements(lon, lat, z = z, radius=100, number=40, time=time_start)
o.seed_elements(lon_ss, lat_ss, z = z_ss, radius=100, number=40, time=time_start)
o.seed_elements(lon_ig, lat_ig, z = z_ig, radius=100, number=40, time=time_start)
o.seed_elements(lon_cn, lat_cn, z = z_cn, radius=100, number=40, time=time_start)
#o.seed_elements(lon_1, lat_1, z = z_1, radius=100, number=40, time=time_start)
o.seed_elements(lon_2, lat_2, z = z_2, radius=100, number=40, time=time_start)
#o.seed_elements(lon_3, lat_3, z = z_3, radius=100, number=40, time=time_start)
o.seed_elements(lon_5, lat_5, z = z_5, radius=100, number=40, time=time_start)

o.set_config('drift:current_uncertainty', .15) 





# Running model (until end of driver data)
#o.run(steps=282, time_step=900, time_step_output=900, outfile='Dev261_9_noangle_nowind_Opendrift.nc') 
o.run(steps=400, time_step=10800)   
#o.run(duration=timedelta(hours=24))
#o.run(steps = 235, time_step = 1800)

# Print and plot results
#print(o)
o.plot(linecolor='z', cmap = 'brg', title = 'Simultaneous seeds at t = 0 for 13/01/2014 to 03/03/2014. Time_step = 3 h')  # Color lines according to depth
o.plot_property('z')
o.animation(filename ='Summer2014.gif')

print(o)
