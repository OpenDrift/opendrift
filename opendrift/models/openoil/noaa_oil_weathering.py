# Oil weathering interface for OpenDrift/OpenOil
# towards NOAA Oil library:
#   https://github.com/NOAA-ORR-ERD/OilLibrary
# Methods below are adapted from PyGnome:
#   https://github.com/NOAA-ORR-ERD/PyGnome

import numpy as np


def mass_transport_coeff(wind_speed):
    c_evap = 0.0025
    mass_transport_coeff = c_evap*np.power(wind_speed, 0.78)
    mass_transport_coeff[wind_speed >= 10] = \
        0.06*c_evap*np.power(wind_speed[wind_speed >= 10], 2)
    return mass_transport_coeff


def evap_decay_constant(substance, wind_speed, sea_water_temperature,
                        area, mass_components):
    K = mass_transport_coeff(wind_speed)  # per element
    f_diff = 1.0
    # vp per element, subcomponent
    vp = np.array([substance.vapor_pressure(t)
                   for t in sea_water_temperature])
    # evaporation expects mw in kg/mol, database is in g/mol
    mw = substance.molecular_weight/1000.
    # sum of mass components, per element
    sum_mi_mw = (mass_components/mw).sum(axis=1)

    gas_constant = 8.314
    decay = (-(area*f_diff*K) / (gas_constant*sea_water_temperature*
                                 sum_mi_mw)).reshape(-1, 1) * vp
    return decay


def water_uptake_coefficient(substance, wind_speed):
    # water uptake rate constant - from database
    K0Y = substance.k0y
    drop_max = 1.0e-5
    k_emul = 6.0 * K0Y * wind_speed * wind_speed / drop_max
    return k_emul
