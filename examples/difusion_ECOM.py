#!/usr/bin/env python
"""
Diffusion_ECOM_EC1 (Verão 2014)
=============
"""
from datetime import datetime, timedelta
from opendrift.readers import reader_ECOM_building_sup
from opendrift.models.oceandrift import OceanDrift
#from opendrift.readers import reader_netCDF_CF_generic



o = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
reader_pcse = reader_ECOM_building_sup.Reader('/home/arian/IC/arquivos_cdf/six_days.nc')
#reader_pcse = reader_ECOM_building.Reader('/home/arian/IC/arquivos_cdf/EC1.cdf')
#reader_pcse = reader_ECOM_building.Reader('/home/arian/IC/arquivos_cdf/z_six_days.nc')


#print(reader_pcse)


o.add_reader([reader_pcse])


time = reader_pcse.start_time

o.seed_elements(lon=-46,lat =-23.9,  number = 200, time=time)
#o.set_config('drift:current_uncertainty', 0)  # 0 is default
o.run(duration=timedelta(hours=130))
#o.run(steps=600, time_step=135)

#o.plot(fast=True)
#o.animation()


o2 = OceanDrift(loglevel=0)  # Set loglevel to 0 for debug information
o2.add_reader([reader_pcse])
o2.seed_elements(lon=-46,lat =-23.9,  number = 200, time=time)
o2.set_config('drift:current_uncertainty', .1) # Difference from first run
o2.run(duration=timedelta(hours=130))
#o2.run(steps=600, time_step=600)
#%%
# Comparing
o2.plot(fast=True)
o2.animation(compare=o, legend=['0.1 m/s std for current components', 'No diffusion'], legend_loc='upper center', filename ='Diffusion_ECOM_2006.gif')
