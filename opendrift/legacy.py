# This is a collection of methods that have been removed from various places in the OpenDrift codebase.


# See earlier version of reader_ROMS_native how the GLS parameters were obtained
def gls_tke(windstress, depth, sea_water_density,
            tke, generic_length_scale, gls_parameters=None):
    '''From LADIM model, based on ROMS files.'''

    g = 9.81
    f0 = 0.1  # mean wave frequency
    c_w = 4.0  # wave mixing parameter
    c_i = 0.2  # coefficient for the interior
    if gls_parameters is None:
        # GLS parameters from ROMS, k-omega closure (see ocean.in)
        p = 0.0
        m = 1.0
        n = 1.0
        cmu0 = 0.5477  # for KANTHA_CLAYSON stability function
    else:
        p = gls_parameters['gls_p']
        m = gls_parameters['gls_m']
        n = gls_parameters['gls_n']
        cmu0 = gls_parameters['gls_cmu0']

    phi = 100. * (windstress/sea_water_density)**(3./2.)

    # dissipation and turbulent length scale for interiour of mixed layer
    eps = cmu0**(3.+p/n)*tke**(3./2.+m/n)*generic_length_scale**(-1./n)
    l_i = c_i * tke**(3./2.) * eps**(-1.)

    # diffusivity for interior of mixed layer
    # c_i = sqrt(2.) * cmu0**3
    ki = c_i * (2.*tke)**0.5 * l_i

    # length scale and diffusivity of wave-enhanced layer
    l_w = np.sqrt(phi / (g*f0))
    kwave = c_w * (2*tke)**0.5 * l_w
    kmix = ki + kwave

    K, N = np.meshgrid(kmix, depths)

    return K

