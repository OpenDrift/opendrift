import numpy as np
import logging


def windspeed_Sundby1983(s):
    windspeed_squared = s.environment.x_wind**2 + s.environment.y_wind**2
    K = 76.1e-4 + 2.26e-4 * windspeed_squared
    #valid = windspeed_squared < 13**2
    return K


def stepfunction(s):
    ''' eddy diffusivity with discontinuity for testing of mixing scheme'''
    K = s.elements.z * 0 + 0.1
    K[s.elements.z < -20] = K[s.elements.z < -20] * 0 + 0.02
    return K
