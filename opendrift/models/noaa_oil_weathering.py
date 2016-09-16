# Oil weathering interface for OpenDrift/OpenOil
# towards NOAA Oil library, adapted from PyGnome

import numpy as np

def mass_transport_coeff(wind_speed):
    c_evap = 0.0025
    mass_transport_coeff = c_evap*np.power(wind_speed, 0.78)
    mass_transport_coeff[wind_speed >= 10] = \
            0.06*c_evap*np.power(wind_speed[wind_speed >= 10], 2)
    return mass_transport_coeff

def evap_decay_constant(substance, wind_speed, sea_water_temperature):
    K = mass_transport_coeff(wind_speed)
    f_diff = 1.0
    vp = np.array([substance.vapor_pressure(t) for t in sea_water_temperature])
    print vp
    print vp.shape
    import sys; sys.exit('noaa')
    
