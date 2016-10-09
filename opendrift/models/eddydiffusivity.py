import numpy as np
import scipy as sp
import logging


def windspeed_Sundby1983(s):
    depths = s.environment_profiles['z']
    windspeed_squared = s.environment.x_wind**2 + s.environment.y_wind**2
    if np.max(windspeed_squared) == 0:
        # Wind not available, check if we can invert from wind stress
        if hasattr(s.environment, 'surface_downward_x_stress'):
            windspeed_squared = s.windspeed_from_stress()**2
    K = 76.1e-4 + 2.26e-4 * windspeed_squared
    # valid = windspeed_squared < 13**2
    Kprofiles, N = sp.meshgrid(K, depths)
    return Kprofiles


def gls_tke(s):
    '''From LADIM model.'''

    if not hasattr(s, 'gls_parameters'):
        logging.info('Searching readers for GLS parameters...')
        for reader_name, reader in s.readers.iteritems():
            if hasattr(reader, 'gls_parameters'):
                s.gls_parameters = reader.gls_parameters
                logging.info('Found gls-parameters in ' + reader_name)
                break  # Success
        if not hasattr(s, 'gls_parameters'):
            logging.info('Did not find gls-parameters in any readers.')
            s.gls_parameters = None

    g = 9.81
    f0 = 0.1  # mean wave frequency
    c_w = 4.0  # wave mixing parameter
    c_i = 0.2  # coefficient for the interior
    if hasattr(s, 'gls_parameters'):
        p = s.gls_parameters['gls_p']
        m = s.gls_parameters['gls_m']
        n = s.gls_parameters['gls_n']
        cmu0 = s.gls_parameters['gls_cmu0']
    else:
        # GLS parameters from ROMS, k-omega closure (see ocean.in)
        p = 0.0
        m = 1.0
        n = 1.0
        cmu0 = 0.5477  # for KANTHA_CLAYSON stability function

    stress = np.sqrt(s.environment.surface_downward_x_stress**2 +
                     s.environment.surface_downward_y_stress**2)

    phi = 100. * (stress/s.sea_water_density())**(3./2.)

    # dissipation and turbulent length scale for interiour of mixed layer
    eps = cmu0**(3.+p/n) * s.environment.turbulent_kinetic_energy**(
        3./2.+m/n) * s.environment.turbulent_generic_length_scale**(-1./n)
    l_i = c_i * s.environment.turbulent_kinetic_energy**(3./2.) * eps**(-1.)

    # diffusivity for interior of mixed layer
    # c_i = sqrt(2.) * cmu0**3
    ki = c_i * (2.*s.environment.turbulent_kinetic_energy)**0.5 * l_i

    # length scale and diffusivity of wave-enhanced layer
    l_w = np.sqrt(phi / (g*f0))
    kwave = c_w * (2*s.environment.turbulent_kinetic_energy)**0.5 * l_w
    kmix = ki + kwave
    depths = s.environment_profiles['z']
    Kprofiles, N = sp.meshgrid(kmix, depths)

    return Kprofiles


def stepfunction(s):
    ''' eddy diffusivity with discontinuity for testing of mixing scheme'''
    Nparticles = range(s.elements.z.shape[0])
    K = s.environment_profiles['z'] * 0 + 0.1
    K[s.environment_profiles['z'] < -20] = \
        K[s.environment_profiles['z'] < -20] * 0 + 0.02
    N, Kprofiles = sp.meshgrid(Nparticles, K)
    return Kprofiles
