import numpy as np
import scipy as sp
import logging



def windspeed_Sundby1983(s):
    depths = s.environment_profiles['z']
    windspeed_squared = s.environment.x_wind**2 + s.environment.y_wind**2
    K = 76.1e-4 + 2.26e-4 * windspeed_squared
    #valid = windspeed_squared < 13**2
    Kprofiles, N = sp.meshgrid(K, depths)
    return Kprofiles


def stepfunction(s):
    ''' eddy diffusivity with discontinuity for testing of mixing scheme'''
    Nparticles = range(s.elements.z.shape[0])
    K = s.environment_profiles['z'] * 0 + 0.1
    K[s.environment_profiles['z'] < -20] = K[s.environment_profiles['z'] < -20] * 0 + 0.02
    N, Kprofiles = sp.meshgrid(Nparticles, K)
    return Kprofiles
