# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2016, Knut-Frode Dagestad, MET Norway

import logging; logger = logging.getLogger(__name__)
from datetime import timedelta
import numpy as np
import scipy as sp
from math import sqrt
import matplotlib.pyplot as plt
import pyproj
import cmocean


def wind_drift_factor_from_trajectory(trajectory_dict, min_period=None):
    '''Estimate wind_drift_fator based on wind and current along given trajectory

    trajectory_dict: dictionary with arrays of same length of the following variables:
    time, lon, lat, x_wind, y_wind, x_sea_water_velocity, y_sea_water_velocity

    Returns array of same length minus one of the fitted wind_drift_factor
    '''

    geod = pyproj.Geod(ellps='WGS84')
    time = trajectory_dict['time']
    try:
        import pandas as pd
        time = pd.to_datetime(time)
    except:
        pass
    if min_period is None:
        ind = range(len(time))
    else:
        timestep = time[1] - time[0]
        s = np.round(min_period.total_seconds()/(timestep).total_seconds()).astype(int)
        ind = np.arange(0, len(time), s).astype(np.int32)
        print('Original timestep (%s) multiplied by %i: %s' % (timestep, s, timestep*s))
        ind2 = ind.copy()
        for i in range(1, s):
            ind2 = np.concatenate((ind2, ind+i))
        ind = ind2
        ind = ind[ind<len(time)]
        time = time[ind]
    cu = trajectory_dict['x_sea_water_velocity'][ind]
    cv = trajectory_dict['y_sea_water_velocity'][ind]
    wu = trajectory_dict['x_wind'][ind]
    wv = trajectory_dict['y_wind'][ind]
    lon = trajectory_dict['lon'][ind]
    lat = trajectory_dict['lat'][ind]
    current_speed = np.sqrt(cu**2+cv**2)
    current_azimuth = np.degrees(np.arctan2(cu, cv))
    wind_speed = np.sqrt(wu**2+wv**2)
    time_step = (time[1]-time[0]).total_seconds()
    # Move forward from each position with current
    lonf, latf, back_az = geod.fwd(lon[0:-1], lat[0:-1], current_azimuth[0:-1], current_speed[0:-1]*time_step)
    # Find distance and azimuth to next observational position
    azimuth_forward, a2, distance = geod.inv(lonf, latf, lon[1:], lat[1:])
    # Find the wind_drift_factor needed to give the remaining distance
    wind_drift_factor = distance / (time_step*wind_speed[0:-1])

    wind_azimuth = np.degrees(np.arctan2(wu, wv))
    azimuth_offset = azimuth_forward - wind_azimuth[0:-1]  # 0 if downwind drift, positive if rightwards drift is needed towards end position
    azimuth_offset = azimuth_offset % 360  # Make sure azimuth angle is between -180 and 180
    azimuth_offset = (azimuth_offset + 360) % 360
    azimuth_offset[azimuth_offset>180] -= 360

    return wind_drift_factor, azimuth_offset

def distance_between_trajectories(lon1, lat1, lon2, lat2):
    '''Calculate the distances [m] between two trajectories'''

    assert len(lon1) == len(lat1) == len(lat1) == len(lat2)
    geod = pyproj.Geod(ellps='WGS84')
    azimuth_forward, a2, distance = geod.inv(lon1, lat1, lon2, lat2)

    return distance

def distance_along_trajectory(lon, lat):
    '''Calculate the distances [m] between points along a trajectory
    assuming the last dimension is the time dimension for the
    trajectory.
    '''

    geod = pyproj.Geod(ellps='WGS84')
    azimuth_forward, a2, distance = geod.inv(lon[...,1:], lat[...,1:], lon[...,:-1], lat[...,:-1])

    return distance

def skillscore_darpa(lon1, lat1, lon2, lat2):
    '''Scoring algorithm used for DARPA float challenge 2021.

    Copied from implementation made by Jean Rabault.
    Assuming 6 positions, separated by 2 days, where first pos (start) are identical.
    '''

    assert len(lon1) == 6
    distances = distance_between_trajectories(lon1, lat1, lon2, lat2)/1000
    if distances[0] != 0:
        logger.warning('DARPA score: first position is not identical '
                       '(distance: %s km)' % distances[0])
    distances = distances[1:]  # For the remaining, we ignore starting position

    dict_DARPA_points = {
        4: 5,
        8: 2,
        16: 1,
        32: 0.25 }
    distance_thresholds = sorted(list(dict_DARPA_points.keys()))

    dict_DARPA_scoring_multiplicator = {
        0: 1,
        1: 2,
        2: 5,
        3: 10,
        4: 20 }

    DARPA_points = 0
    for crrt_ind, crrt_distance in enumerate(distances):
        if crrt_distance >= distance_thresholds[-1]:
            break

        for crrt_threshold in distance_thresholds:
            if crrt_distance < crrt_threshold:
                crrt_points = dict_DARPA_points[crrt_threshold]
                break

        crrt_multiplicator = dict_DARPA_scoring_multiplicator[crrt_ind]

        DARPA_points += int(crrt_points * crrt_multiplicator)

    return DARPA_points

def skillscore_liu_weissberg(lon_obs, lat_obs, lon_model, lat_model, tolerance_threshold=1):
    ''' calculate skill score from normalized cumulative seperation
    distance. Liu and Weisberg 2011.

    Returns the skill score bewteen 0. and 1.
    '''

    lon_obs = np.array(lon_obs)
    lat_obs = np.array(lat_obs)
    lon_model = np.array(lon_model)
    lat_model = np.array(lat_model)
    d = distance_between_trajectories(lon_obs, lat_obs, lon_model, lat_model)
    l = distance_along_trajectory(lon_obs, lat_obs)
    s = np.divide(np.nansum(d, axis = -1), np.nansum(np.nancumsum(l, axis = -1), axis = -1))
    if tolerance_threshold==0:
        skillscore = 0
    else:
        skillscore = np.maximum(0, 1 - s/tolerance_threshold)

    return skillscore

def plot_wind_drift_factor(wdf, azimuth, wmax_plot=None):
    '''Polar plot of array of wind drift factor, with associated azimuthal offset'''

    wmax = wdf.max()
    wbins = np.arange(0, wmax+.005, .005)
    abins = np.linspace(-180, 180, 30)
    hist, _, _ = np.histogram2d(azimuth, wdf, bins=(abins, wbins))
    A, W = np.meshgrid(abins, wbins)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location('N', offset=0)
    ax.set_theta_direction(-1)
    pc = ax.pcolormesh(np.radians(A), W, hist.T, cmap=cmocean.cm.dense)
    plt.arrow(np.pi, wmax, 0, -wmax, width=.015, facecolor='k', zorder=100,
              head_width=.8, lw=2, head_length=.005, length_includes_head=True)
    plt.text(np.pi, wmax*.5, ' Wind direction', fontsize=18)
    plt.grid()
    if wmax_plot is not None:
        plt.ylim([0, wmax_plot])
    plt.show()


def oil_wave_entrainment_rate_li2017(dynamic_viscosity, oil_density, interfacial_tension,
                                     significant_wave_height=None, wave_breaking_fraction=None,
                                     wind_speed=None, sea_water_density=1028.):
    # Z. Li, M.L. Spaulding, D. French McCay, J. Mar. Pollut. Bull. (2016):
    # An algorithm for modeling entrainment and naturally and chemically dispersed
    # oil droplet size distribution under surface breaking wave conditions

    if wave_breaking_fraction is None:
        if wind_speed is None:
            raise ValueError('wave_breaking_fraction or wind_speed must be provided')
        wave_breaking_fraction = wave_breaking_fraction_from_wind(wind_speed)
    if significant_wave_height is None:
        if wind_speed is None:
            raise ValueError('significant_wave_height or wind_speed must be provided')
        significant_wave_height = significant_wave_height_from_wind_neumann_pierson(wind_speed)
    g = 9.81
    delta_rho = sea_water_density - oil_density
    d_o = 4*np.sqrt(interfacial_tension / (delta_rho*g))
    we = sea_water_density*g*significant_wave_height*d_o/interfacial_tension
    oh = dynamic_viscosity/np.sqrt(oil_density*interfacial_tension*d_o)
    with np.errstate(divide='ignore'):
        entrainment_rate = (4.604e-10*we**1.805*oh**-1.023)*wave_breaking_fraction
    return entrainment_rate

def significant_wave_height_from_wind_neumann_pierson(wind_speed):
    # Neumann and Pierson, 1966
    # WMO 1998
    return 0.0246*np.power(wind_speed, 2)

def wave_breaking_fraction_from_wind(wind_speed, wave_period=None):
    # TODO: We should also have an option here for
    # the case when wave height is given, but no wind
    if wave_period is None:
        wave_period = wave_period_from_wind(wind_speed)
    f = 0.032*(wind_speed - 5)/wave_period
    f[f < 0] = 0
    return f

def wave_period_from_wind(wind_speed):
    # Pierson-Moskowitz if period not available from readers
    # WMO guide to wave analysis and forecasting pp. 14, WMO (1998)
    # fallback value if no wind speed or Hs to avoid division by zero
    wind_speed = np.atleast_1d(wind_speed)
    omega = 5*np.ones(wind_speed.shape)  # Angular frequency
    omega[wind_speed>0] = 0.877*9.81/(1.17*wind_speed[wind_speed>0])
    return 2*np.pi/omega

def verticaldiffusivity_Sundby1983(windspeed, depth, mixedlayerdepth=50, background_diffusivity=0):
    ''' Vertical diffusivity from Sundby (1983)

    S. Sundby (1983): A one-dimensional model for the vertical
        distribution of pelagic fish eggs in the mixed layer
        Deep Sea Research (30) pp. 645-661
    '''

    K = 76.1e-4 + 2.26e-4 * windspeed*windspeed * np.ones(np.atleast_1d(depth.shape))
    K[depth>mixedlayerdepth-1] = (K[depth>mixedlayerdepth-1]+background_diffusivity)/2  # Transition
    K[depth>=mixedlayerdepth] = background_diffusivity  # Cutoff below mixed layer
    # valid = windspeed_squared < 13**2
    return K

def verticaldiffusivity_Large1994(windspeed, depth, mixedlayerdepth=50, background_diffusivity=0):
    ''' Vertical diffusivity from Large et al. (1994)

    Depending on windspeed, depth and mixed layer depth (default 50m).'''

    # Defining two helper methods:

    def stabilityfunction(sigma):
        # controls stratification regimes for diffusivity
        # a value of 1 represent neutrally buoyant conditions (no stratification)
        # a value between 0 and 1 represents stable stratification
        return 0.2 # should be depth dependent, too

    def G(sigma):
        # vertical shape function for eddy diffusivity
        a1 = 1.
        a2 = -2
        a3 = 1
        G = a1*sigma + a2*sigma**2 + a3*sigma**3
        below = np.where(G>=1)
        G[below]=G[below]*0.
        return G

    depth = np.abs(depth)  # Making sure depth is positive value
    MLD = mixedlayerdepth  # shorthand
    rhoa = 1.22  # Air density
    cd = 1.25e-3  # Kara et al. 2007
    windstress = windspeed*windspeed * cd * rhoa

    K = MLD * stabilityfunction(depth/MLD) * 0.4 * G(depth/MLD) * windstress + (depth/MLD)*background_diffusivity
    K[depth>=MLD] = background_diffusivity

    return K

def verticaldiffusivity_stepfunction(depth, MLD=20,
                                     k_above=.1, k_below=.02):
    ''' eddy diffusivity with discontinuity for testing of mixing scheme'''
    depth = np.abs(depth)  # Making sure depth is positive value
    K = k_above*np.ones(depth.shape)
    K[depth>MLD] = k_below
    return K

def gls_tke(windstress, depth, sea_water_density,
            tke, generic_length_scale, gls_parameters=None):
    '''From LADIM model.'''

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

def plot_stokes_profile(profiles, view=['vertical', 'birdseye', 'u', 'v'], filename=None):
    '''Plot vertical profile of Stokes drift

    Args:
        list of dictionary with:
            u, v: components of Stokes drift vector
            z: depth in meters, 0 at surface and negative below (e.g. -5 = 5m depth)
            kwargs: forwarded to plot method
        view:
            - 'vertical': magnitude versus depth
            - 'birdseye': viewed from above
            - or list of both above with one axis for each case (default)
    '''

    if len(view)>2:
        fig, axs = plt.subplots(2, 2, figsize=(12,12))
        axs = axs.ravel()
    else:
        fig, axs = plt.subplots(1, len(view))
    for vi, ax in zip(view, axs):
        for p in profiles:
            if 'kwargs' in p:
                kwargs = p['kwargs']
            else:
                kwargs=None
            u = p['u']
            v = p['v']
            z = p['z']
            if vi == 'vertical':
                speed = np.sqrt(u**2 + v**2)
                ax.plot(speed, -z, **kwargs)
            elif vi == 'birdseye':
                ax.plot(u, v, **kwargs)
            elif vi == 'u':
                ax.plot(u, -z, **kwargs)
            elif vi == 'v':
                ax.plot(v, -z, **kwargs)

        if vi == 'vertical':
            ax.invert_yaxis()
            ax.set_xlabel('Speed  [m/s]')
            ax.set_ylabel('Depth  [m]')
            ax.set_ylim(ax.get_ylim()[0], 0)
        elif vi == 'birdseye':
            ax.set_xlabel('u  [m/s]')
            ax.set_ylabel('v  [m/s]')
            ax.set_aspect('equal')
            xlim = np.array(ax.get_xlim())
            ylim = np.array(ax.get_ylim())
            m = np.maximum(np.abs(xlim).max(), np.abs(ylim).max())
            ax.set_xlim(-m, m)
            ax.set_ylim(-m, m)
        if vi == 'u':
            ax.invert_yaxis()
            ax.set_xlabel('U-component  [m/s]')
            ax.set_ylabel('Depth  [m]')
            ax.set_ylim(ax.get_ylim()[0], 0)
        if vi == 'v':
            ax.invert_yaxis()
            ax.set_xlabel('V-component  [m/s]')
            ax.set_ylabel('Depth  [m]')
            ax.set_ylim(ax.get_ylim()[0], 0)
        ax.legend()

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)

def stokes_transport_monochromatic(mean_wave_period, significant_wave_height):
    mean_wave_frequency = 2.*np.pi/mean_wave_period
    return mean_wave_frequency * np.power(significant_wave_height, 2) / 16

def stokes_drift_profile_monochromatic(stokes_u_surface, stokes_v_surface,
                                       significant_wave_height, mean_wave_period, z):
    """
    Vertical Stokes drift profile assuming a single monochromatic wave
    Breivik, O., Janssen, P., Bidlot, J., 2014. Approximate stokes drift profiles in deep water.
    J. Phys. Oceanogr. 44, 2433-2445.  doi:10.1175/JPO-D-14-0020.1.
    """
    stokes_surface_speed = np.sqrt(stokes_u_surface**2 +
                                   stokes_v_surface**2)

    km = stokes_surface_speed / (
            2*stokes_transport_monochromatic(mean_wave_period, significant_wave_height))

    stokes_speed = stokes_surface_speed*np.exp(2*km*z)

    zeromask = stokes_surface_speed == 0
    stokes_u = stokes_speed*stokes_u_surface/stokes_surface_speed
    stokes_v = stokes_speed*stokes_v_surface/stokes_surface_speed
    stokes_u[zeromask] = 0
    stokes_v[zeromask] = 0

    return stokes_u, stokes_v, stokes_speed

def stokes_drift_profile_exponential(stokes_u_surface, stokes_v_surface,
                                     significant_wave_height, mean_wave_period, z):

    """
    Breivik, O., Janssen, P., Bidlot, J., 2014. Approximate stokes drift profiles in deep water.
    J. Phys. Oceanogr. 44, 2433-2445.  doi:10.1175/JPO-D-14-0020.1.
    """

    stokes_surface_speed = np.sqrt(stokes_u_surface**2 +
                                     stokes_v_surface**2)

    km = stokes_surface_speed / (
          2*stokes_transport_monochromatic(mean_wave_period, significant_wave_height))
    ke = km/3

    stokes_speed = stokes_surface_speed*np.exp(2*ke*z)/(1-8*ke*z)

    zeromask = stokes_surface_speed == 0
    stokes_u = stokes_speed*stokes_u_surface/stokes_surface_speed
    stokes_v = stokes_speed*stokes_v_surface/stokes_surface_speed
    stokes_u[zeromask] = 0
    stokes_v[zeromask] = 0

    return stokes_u, stokes_v, stokes_speed

def stokes_drift_profile_phillips(stokes_u_surface, stokes_v_surface,
                                  significant_wave_height, mean_wave_period, z):
    """
    Calculate vertical Stokes drift profile from
    Breivik et al. 2016, A Stokes drift approximation based on the Phillips spectrum, Ocean Mod. 100
    """

    alpha = 0.0083
    stokes_surface_speed = np.sqrt(stokes_u_surface**2 +
                                   stokes_v_surface**2)

    beta = 1
    km = stokes_surface_speed * (1-2*beta/3)/ (
            2*stokes_transport_monochromatic(mean_wave_period, significant_wave_height))

    stokes_speed = stokes_surface_speed*(np.exp(2*km*z) -
        beta*np.sqrt(2*np.pi*km*np.abs(z))*sp.special.erfc(np.sqrt(2*km*np.abs(z))))

    zeromask = stokes_surface_speed == 0
    stokes_u = stokes_speed*stokes_u_surface/stokes_surface_speed
    stokes_v = stokes_speed*stokes_v_surface/stokes_surface_speed
    stokes_u[zeromask] = 0
    stokes_v[zeromask] = 0

    return stokes_u, stokes_v, stokes_speed

def stokes_drift_profile_windsea_swell(stokes_u_surface, stokes_v_surface,
                                       swell_mean_direction_to, swell_mean_period, swell_height,
                                       wind_sea_mean_direction_to, wind_sea_mean_period, wind_sea_height, z):
    """
    Calculate vertical Stokes drift profile from
    Breivik, O., and K. H. Christensen, 2020: A Combined Stokes Drift Profile under Swell and Wind Sea.
        J. Phys. Oceanogr., 50, 2819-2833, https://doi.org/10.1175/JPO-D-20-0087.1.
    """

    # NB / TODO: assuming here that u is east and v is north

    # Calculate swell Stokes component
    th_ws_N = np.cos(np.radians(wind_sea_mean_direction_to))  # unit vector
    th_ws_E = np.sin(np.radians(wind_sea_mean_direction_to))  # unit vector
    th_sw_N = np.cos(np.radians(swell_mean_direction_to))     # unit vector
    th_sw_E = np.sin(np.radians(swell_mean_direction_to))     # unit vector
    stokes_swell_surface_speed = (stokes_u_surface*th_ws_N - stokes_v_surface*th_ws_E) / \
        (th_sw_E*th_ws_N - th_sw_N*th_ws_E)
    stokes_swell_surface_u = stokes_swell_surface_speed*th_sw_E
    stokes_swell_surface_v = stokes_swell_surface_speed*th_sw_N

    stokes_swell_u, stokes_swell_v, stokes_swell_speed = stokes_drift_profile_monochromatic(
        stokes_swell_surface_u, stokes_swell_surface_v,
        swell_height, swell_mean_period, z)

    # Calculate wind Stokes component, with constraint that
    # wind and swell Stokes at surface sums to total
    stokes_wind_surface_u = stokes_u_surface - stokes_swell_surface_u
    stokes_wind_surface_v = stokes_v_surface - stokes_swell_surface_v

    stokes_wind_u, stokes_wind_v, stokes_wind_speed = stokes_drift_profile_phillips(
        stokes_wind_surface_u, stokes_wind_surface_v, wind_sea_height, wind_sea_mean_period, z)

    stokes_u = stokes_swell_u + stokes_wind_u
    stokes_v = stokes_swell_v + stokes_wind_v
    stokes_speed = np.sqrt(stokes_u**2+stokes_v**2)

    return stokes_u, stokes_v, stokes_speed


def ftle(X, Y, delta, duration):
    """Calculate Finite Time Lyapunov Exponents"""
    # From Johannes Rohrs
    nx = X.shape[0]
    ny = X.shape[1]
    J = np.empty([nx,ny,2,2], np.float32)
    FTLE = np.empty([nx,ny], np.float32)

    # gradient
    dx = np.gradient(X)
    dy = np.gradient(Y)

    # Jacobian
    J[:,:,0,0] = dx[0] / (2*delta)
    J[:,:,1,0] = dy[0] / (2*delta)
    J[:,:,0,1] = dx[1] / (2*delta)
    J[:,:,1,1] = dy[1] / (2*delta)

    for i in range(0,nx):
        for j in range(0,ny):
            # Green-Cauchy tensor
            D = np.dot(np.transpose(J[i,j]), J[i,j])
            # its largest eigenvalue
            lamda = np.linalg.eigvals(D)
            FTLE[i,j] = np.log(np.sqrt(max(lamda)))/np.abs(duration)

    return FTLE

__stokes_coefficients__ = None

def wave_stokes_drift_parameterised(wind, fetch):
    """
    Parameterise stokes drift based on pre calculated tables and fetch.
    """
    global __stokes_coefficients__

    Wf       = {}
    Wf['5000'] = (0.0173,0.0160,0.0152,0.0145,0.0139,0.0135,
                0.0132,0.0129,0.0126,0.0124,0.0122,0.0121,
                0.0119,0.0118,0.0117,0.0116,0.0114,0.0113,
                0.0112,0.0112,0.0111,0.0110,0.0109,0.0109,
                0.0108,0.0107,0.0106,0.0106,0.0106,0.0105)

    Wf['25000'] = (0.0173,0.0197,0.0201,0.0185,0.0181,0.0176,
                    0.0171,0.0167,0.0164,0.0160,0.0158,0.0155,
                    0.0153,0.0151,0.0149,0.0147,0.0146,0.0144,
                    0.0143,0.0142,0.0140,0.0139,0.0138,0.0137,
                    0.0136,0.0135,0.0135,0.0134,0.0133,0.0132)

    Wf['50000'] = (0.0173,0.0197,0.0210,0.0216,0.0201,0.0194,
                    0.0190,0.0186,0.0183,0.0179,0.0176,0.0173,
                    0.0171,0.0168,0.0166,0.0164,0.0162,0.0160,
                    0.0159,0.0157,0.0156,0.0155,0.0153,0.0152,
                    0.0151,0.0150,0.0149,0.0148,0.0147,0.0146)
    n = {'5000': 3, '25000':6, '50000': 6}  # Polynom order

    if __stokes_coefficients__ is None:
        c = {}
        for f in ['5000', '25000', '50000']:
            c[f] = np.polyfit(
                range(len(Wf[f])), Wf[f], n[f])
        __stokes_coefficients__ = c

    windspeed = np.sqrt(wind[0]**2 + wind[1]**2)
    windspeed[windspeed>30] = 30
    wf = np.polyval(__stokes_coefficients__[fetch], windspeed)
    stokes_drift_x_velocity = wind[0]*wf
    stokes_drift_y_velocity = wind[1]*wf

    return stokes_drift_x_velocity, stokes_drift_y_velocity

__stokes_hs_coefficients__ = None
def wave_significant_height_parameterised(wind, fetch):
    """
    Parameterise significant wave height based on pre calculated tables and fetch.
    """

    global __stokes_hs_coefficients__

    Sw       = {}
    Sw['5000'] = (0.030,0.077,0.124,0.170,0.216,0.263,
                0.311,0.360,0.409,0.459,0.509,0.560,
                0.612,0.664,0.716,0.771,0.823,0.876,
                0.932,0.987,1.041,1.095,1.152,1.210,
                1.265,1.319,1.375,1.434,1.494,1.552)

    Sw['25000'] = (0.030,0.122,0.251,0.336,0.442,0.546,
                    0.650,0.753,0.856,0.959,1.063,1.168,
                    1.273,1.379,1.486,1.593,1.702,1.811,
                    1.920,2.030,2.142,2.254,2.366,2.478,
                    2.592,2.707,2.822,2.936,3.051,3.166)

    Sw['50000'] = (0.030,0.122,0.274,0.474,0.591,0.724,
                    0.873,1.021,1.168,1.314,1.460,1.606,
                    1.752,1.898,2.045,2.192,2.340,2.489,
                    2.639,2.789,2.940,3.092,3.244,3.397,
                    3.551,3.706,3.862,4.017,4.173,4.330)

    if __stokes_hs_coefficients__ is None:
        c = {}
        for f in ['5000', '25000', '50000']:
            c[f] = np.polyfit(
                range(len(Sw[f])), Sw[f], 1)
        __stokes_hs_coefficients__ = c

    windspeed = np.sqrt(wind[0]**2 + wind[1]**2)
    windspeed[windspeed>30] = 30
    wave_significant_height = np.polyval(
        __stokes_hs_coefficients__[fetch], windspeed)

    return wave_significant_height

class PhysicsMethods:
    """Physics methods to be inherited by OpenDriftSimulation class"""

    @staticmethod
    def sea_water_density(T=10., S=35.):
        '''The function gives the density of seawater at one atmosphere
        pressure as given in :

        N.P. Fofonoff and R.C. Millard Jr.,1983,
        Unesco technical papers in marine science no. 44.

        S   = Salinity in promille of the seawater
        T   = Temperature of the seawater in degrees Celsius
        '''

        if np.atleast_1d(T).max() > 100:
            raise ValueError('Temperature should be in celcius, but is > 100')

        R4 = 4.8314E-04
        DR350 = 28.106331

        # Pure water density at atmospheric pressure
        # Bigg P.H. (1967) BR. J. Applied Physics pp.:521-537

        R1 = ((((6.536332E-09 * T - 1.120083E-06) * T + 1.001685E-04) *
              T - 9.095290E-03) * T + 6.793952E-02) * T - 28.263737

        # Seawater density at atmospheric pressure
        # coefficients involving salinity :

        R2 = (((5.3875E-09 * T - 8.2467E-07) * T + 7.6438E-05) *
              T - 4.0899E-03) * T + 8.24493E-01

        R3 = (-1.6546E-06*T+1.0227E-04)*T-5.72466E-03

        # International one-atmosphere equation of state of seawater :

        SIG = R1 + (R4*S + R3*np.sqrt(S) + R2)*S
        Dens0 = SIG + DR350 + 1000.
        return Dens0


    def skillscore_trajectory(self, lon_obs, lat_obs, time_obs, duration=None,
                              max_time_offset=timedelta(seconds=60), method='liu-weissberg',
                              **kwargs):
        '''Calculate skill-score comparing simulated and observed trajectories

        A skill score is calculated for each of the modelled trajectories,
        starting from the closest time/obs of the given modelled trajectory

        Parameters
        ----------
        lon_obs : array_like
        lat_obs : array_like
            The longitude and latitudes of the observed trajectory
        time_obs : array of datetime
            The time corresponding to the observed positions
        duration : timedelta
            If None, skill score is calculated for the full trajectory, from the start
            Otherwise, only this sub-part of the trajectory is used.
        max_time_offset : timedelta
            The maximum allowed time difference between observed and modelled time
            Default: 1 minute
        method : string
            Currently available: 'liu-wessberg'
            To be implemented: 'darpa'
        **kwargs : further arguments for skillscore method (e.g. tolerance_threshold)

        Returns
        -------
        skillscore : array_like
            One value for each trajectory
        '''

        time_obs = np.array(time_obs)
        timesteps_obs = time_obs[1:] - time_obs[0:-1]
        if timesteps_obs.min() != timesteps_obs.max():
            raise ValueError('Time step of trajectory is variable: %s to %s' %
                             (timesteps_obs.min(), timesteps_obs.max()))
        if timesteps_obs.min() != self.time_step_output:
            raise ValueError('Time step of trajectory is %s, but time step of simulation is %s.'
                    % (timesteps_obs.min(), self.time_step_output))

        from bisect import bisect_left
        times = np.array(self.get_time_array()[0])
        skillscore = np.zeros(self.num_elements_total())
        for i in range(self.num_elements_total()):
            lon = self.history['lon'][i]
            lat = self.history['lat'][i]
            status = self.history['status'][i]
            lon_model = lon[status==0]
            lat_model = lat[status==0]
            time_model = times[status==0]

            indo = bisect_left(time_obs, time_model[0], hi=len(time_obs)-1)
            if time_obs[indo+1] - time_model[0] < time_obs[indo] - time_model[0]:
                indo = indo + 1
            # Check that modeled and observed (sub)trajectories start at same position
            np.testing.assert_almost_equal(lon_model[0], lon_obs[indo], 3)
            np.testing.assert_almost_equal(lat_model[0], lat_obs[indo], 3)

            time_offset = np.abs((time_model[0] - time_obs[indo]).total_seconds())
            if time_offset > max_time_offset.total_seconds():
                raise ValueError('Too large offset: %s' % time_offset)
            length = len(time_obs) - indo
            indo = slice(indo, indo+length)
            skillscore[i] = skillscore_liu_weissberg(lon_obs[indo], lat_obs[indo],
                    lon_model[0:length], lat_model[0:length], **kwargs)

        return skillscore



    def advect_ocean_current(self, factor=1):

        cdf = self.elements.current_drift_factor
        cdfmin = cdf.min()
        cdfmax = cdf.max()
        if cdfmin != 1 or cdfmax != 1:
            if cdfmin == cdfmax:
                logger.debug('Using currrent drift factor of %s' % cdf)
            else:
                logger.debug('Using currrent drift factor between %s and %s'
                             % (cdfmin, cdfmax))
        factor = factor*cdf
        # Runge-Kutta scheme
        if self.get_config('drift:advection_scheme')[0:11] == 'runge-kutta':
            x_vel = self.environment.x_sea_water_velocity
            y_vel = self.environment.y_sea_water_velocity

            # Find midpoint
            az = np.degrees(np.arctan2(x_vel, y_vel))
            speed = np.sqrt(x_vel*x_vel + y_vel*y_vel)
            dist = speed*self.time_step.total_seconds()*.5
            geod = pyproj.Geod(ellps='WGS84')
            mid_lon, mid_lat, dummy = geod.fwd(self.elements.lon,
                                               self.elements.lat,
                                               az, dist, radians=False)
            # Find current at midpoint, a half timestep later
            logger.debug('Runge-kutta, fetching half time-step later...')
            mid_env, profiles, missing = self.env.get_environment(
                ['x_sea_water_velocity', 'y_sea_water_velocity'],
                self.time + self.time_step/2,
                mid_lon, mid_lat, self.elements.z, profiles=None)
            if self.get_config('drift:advection_scheme') == 'runge-kutta4':
                logger.debug('Runge-kutta 4th order...')
                x_vel2 = mid_env['x_sea_water_velocity']
                y_vel2 = mid_env['y_sea_water_velocity']
                az2 = np.degrees(np.arctan2(x_vel2, y_vel2))
                speed2 = np.sqrt(x_vel2*x_vel2 + y_vel2*y_vel2)
                dist2 = speed2*self.time_step.total_seconds()*.5
                lon2, lat2, dummy = \
                    geod.fwd(self.elements.lon,
                             self.elements.lat,
                             az2, dist2, radians=False)
                env2, profiles, missing = self.env.get_environment(
                    ['x_sea_water_velocity', 'y_sea_water_velocity'],
                    self.time + self.time_step/2,
                    lon2, lat2, self.elements.z, profiles=None)
                # Third step
                x_vel3 = env2['x_sea_water_velocity']
                y_vel3 = env2['y_sea_water_velocity']
                az3 = np.degrees(np.arctan2(x_vel3, y_vel3))
                speed3 = np.sqrt(x_vel3*x_vel3 + y_vel3*y_vel3)
                dist3 = speed3*self.time_step.total_seconds()*.5
                lon3, lat3, dummy = \
                    geod.fwd(self.elements.lon,
                             self.elements.lat,
                             az3, dist3, radians=False)
                env3, profiles, missing = self.env.get_environment(
                    ['x_sea_water_velocity', 'y_sea_water_velocity'],
                    self.time + self.time_step,
                    lon3, lat3, self.elements.z, profiles=None)
                # Fourth step
                x_vel4 = env3['x_sea_water_velocity']
                y_vel4 = env3['y_sea_water_velocity']
                u4 = (x_vel + 2*x_vel2 + 2* x_vel3 + x_vel4)/6.0
                v4 = (y_vel + 2*y_vel2 + 2* y_vel3 + y_vel4)/6.0
                # Move particles using runge-kutta4 velocity
                self.update_positions(u4*factor, v4*factor)

            else:
                # Move particles using runge-kutta velocity
                self.update_positions(
                        factor*mid_env['x_sea_water_velocity'],
                        factor*mid_env['y_sea_water_velocity'])
        elif self.get_config('drift:advection_scheme') == 'euler':
            # Euler scheme
            self.update_positions(
                    factor*self.environment.x_sea_water_velocity,
                    factor*self.environment.y_sea_water_velocity)
        else:
            raise ValueError('Drift scheme not recognised: ' +
                             self.get_config('drift:advection_scheme'))

    def advect_with_sea_ice(self, factor=1):
        if hasattr(self.environment, 'sea_ice_x_velocity'):
            self.update_positions(
                    factor*self.environment.sea_ice_x_velocity,
                    factor*self.environment.sea_ice_y_velocity)
        else:
            if not hasattr(self.environment, 'x_sea_water_velocity'):
                logger.info('No sea ice velocity available')
                return
            # Sea ice velocity, rule-of-thumb from Nordam,
            # doi:10.1016/j.marpolbul.2019.01.019
            ice_velocity_x = self.environment.x_sea_water_velocity + \
                    0.015*self.environment.x_wind
            ice_velocity_y = self.environment.y_sea_water_velocity + \
                    0.015*self.environment.y_wind
            self.update_positions(
                    factor*ice_velocity_x,
                    factor*ice_velocity_y)

    def advect_wind(self, factor=1):
        # Elements at/near ocean surface (z>wind_drift_depth) are advected with given percentage
        # of wind speed. NB: Only basic Euler scheme is implemented

        # Make copy to prevent outside value to be modified
        wind_drift_factor = self.elements.wind_drift_factor.copy()

        # Convert wind_drift_factor to array
        if len(np.atleast_1d(wind_drift_factor)) == 1:
            wind_drift_factor = wind_drift_factor*np.ones(len(self.elements))
        # Convert z to array
        if not hasattr(self.elements, 'z'):
            z = np.zeros(len(self.elements))  # Assumed on surface, if no z
        else:
            z = self.elements.z
        if len(np.atleast_1d(z)) == 1:
            z = z*np.ones(len(self.elements))

        try:
            wind_drift_depth = self.get_config('drift:wind_drift_depth')
        except:
            wind_drift_depth = 0
        if wind_drift_depth == 0:
            surface_only = True
        else:
            wind_drift_depth = np.abs(wind_drift_depth)*np.ones(len(self.elements))
            surface_only = False

        surface = self.elements.z >= -wind_drift_depth
        if surface.sum() == 0:
            if surface_only is True:
                logger.debug('All elements are below surface, no wind-induced shear drift')
            else:
                logger.debug('All elements are below %fm, no wind-induced shear drift' % wind_drift_depth[0])
            return

        wdf = wind_drift_factor.copy()
        wdf_air = wdf.copy()
        if surface_only is False:
            # linear decrease from surface down to wind_drift_depth
            wdf = wdf*(wind_drift_depth+self.elements.z)/wind_drift_depth
            # Resetting wdf for elements in air
            wdf[self.elements.z>0] = wdf_air[self.elements.z>0]
        wdf[~surface] = 0.0
        wdfmin = wdf[surface].min()
        wdfmax = wdf[surface].max()

        x_wind = self.environment.x_wind.copy()
        y_wind = self.environment.y_wind.copy()

        try:
            if self.get_config('drift:relative_wind') is True:
                # Use wind relative to (subtracted) ocean current
                x_wind = x_wind - self.environment.x_sea_water_velocity
                y_wind = y_wind - self.environment.y_sea_water_velocity
        except:
            pass

        speed = np.sqrt(x_wind[surface]*x_wind[surface] + y_wind[surface]*y_wind[surface])
        if wdf[surface].max() == 0:
            logger.debug('Wind drift factor is 0')
            return
        if speed.max() == 0:
            logger.debug('No wind for wind-sheared ocean drift')
            return

        speed = speed*wdf[surface]
        if surface_only is True:
            logger.debug('Advecting %s of %i elements at surface with '
                          'wind-sheared ocean current (%f m/s - %f m/s)'
                          % (np.sum(surface), self.num_elements_active(),
                             speed.min(), speed.max()))
        else:
            logger.debug('Advecting %s of %i elements above %.3fm with '
                          'wind-sheared ocean current (%f m/s - %f m/s)'
                          % (np.sum(surface), self.num_elements_active(),
                             wind_drift_depth[0],
                             speed.min(), speed.max()))

        self.update_positions(x_wind*wdf*factor, y_wind*wdf*factor)

    def stokes_drift(self, factor=1):

        if self.get_config('drift:stokes_drift') is False:
            logger.debug('Stokes drift not activated')
            return

        if np.max(np.array(
            self.environment.sea_surface_wave_stokes_drift_x_velocity+
            self.environment.sea_surface_wave_stokes_drift_y_velocity)) \
                == 0:
            logger.debug('No Stokes drift velocity available')
            return

        wave_height = self.significant_wave_height()
        wave_period = self.wave_period()
        if np.max(np.array(wave_height)) == 0:
            logger.debug('Stokes drift is available, but not Hs: using Hs=1 for Stokes profile')
            wave_height = 1
        if np.max(np.array(wave_period)) == 0:
            logger.debug('Stokes drift is available, but not Tp: using Tp=8 for Stokes profile')
            wave_period = 8

        stokes_profile = self.get_config('drift:stokes_drift_profile')
        if stokes_profile == 'monochromatic':
            stokes_u, stokes_v, s = stokes_drift_profile_monochromatic(
                self.environment.sea_surface_wave_stokes_drift_x_velocity,
                self.environment.sea_surface_wave_stokes_drift_y_velocity,
                wave_height, wave_period, self.elements.z)
        if stokes_profile == 'exponential':
            stokes_u, stokes_v, s = stokes_drift_profile_exponential(
                self.environment.sea_surface_wave_stokes_drift_x_velocity,
                self.environment.sea_surface_wave_stokes_drift_y_velocity,
                wave_height, wave_period, self.elements.z)
        if stokes_profile == 'Phillips':
            stokes_u, stokes_v, s = stokes_drift_profile_phillips(
                self.environment.sea_surface_wave_stokes_drift_x_velocity,
                self.environment.sea_surface_wave_stokes_drift_y_velocity,
                wave_height, wave_period, self.elements.z)
        elif stokes_profile == 'windsea_swell':
            stokes_u, stokes_v, s = stokes_drift_profile_windsea_swell(
                self.environment.sea_surface_wave_stokes_drift_x_velocity,
                self.environment.sea_surface_wave_stokes_drift_y_velocity,
                swell_mean_direction_to = self.environment.sea_surface_swell_wave_to_direction,
                swell_mean_period = self.environment.sea_surface_swell_wave_peak_period_from_variance_spectral_density,
                swell_height = self.environment.sea_surface_swell_wave_significant_height,
                wind_sea_mean_direction_to = self.environment.sea_surface_wind_wave_to_direction,
                wind_sea_mean_period = self.environment.sea_surface_wind_wave_mean_period,
                wind_sea_height = self.environment.sea_surface_wind_wave_significant_height,
                z=self.elements.z)

        self.update_positions(stokes_u*factor, stokes_v*factor)
        if s.min() == s.max():
            logger.debug('Advecting with Stokes drift (%s m/s)' % s.min())
        else:
            logger.debug('Advecting with Stokes drift (%s to %s m/s)' %
                      (s.min(), s.max()))

    def resurface_elements(self, minimum_depth):
        # Keep surfacing elements in water column as default,
        # i.e. no formation of surface slick
        surface = np.where(self.elements.z >= 0)[0]
        self.elements.z[surface] = minimum_depth

    def calculate_missing_environment_variables(self):

        # Missing significant wave height
        if hasattr(self.environment,
                   'sea_surface_wave_significant_height') and \
                self.environment.sea_surface_wave_significant_height.max() == 0:
            Hs = self.significant_wave_height()
            logger.debug('Calculating Hs from wind, min: %f, mean: %f, max: %f' %
                          (Hs.min(), Hs.mean(), Hs.max()))

        # Missing wave periode
        if hasattr(self.environment,
                   'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment') and \
                self.environment.sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment.max() == 0:
            wave_period = self.wave_period()
            logger.debug('Calculating wave period from wind, min: %f, mean: %f, max: %f' %
                          (wave_period.min(), wave_period.mean(), wave_period.max()))


    def wind_speed(self):
        return np.sqrt(self.environment.x_wind**2 +
                       self.environment.y_wind**2)

    def current_speed(self):
        return np.sqrt(self.environment.x_sea_water_velocity**2 +
                       self.environment.y_sea_water_velocity**2)

    def significant_wave_height(self):
        # Significant wave height, parameterise from wind if not available
        if hasattr(self.environment,
                   'sea_surface_wave_significant_height') and \
                self.environment.sea_surface_wave_significant_height.max() > 0:
            Hs = self.environment.sea_surface_wave_significant_height
        else:
            ## Neumann and Pierson, 1966
            #return 0.2*np.power(self.wind_speed(), 2)/9.81
            # WMO 1998
            Hs = 0.0246*np.power(self.wind_speed(), 2)
            self.environment.sea_surface_wave_significant_height = Hs

        return Hs

    def _wave_frequency(self):
        # Note: this is angular frequency, 2*pi*fp
        # Pierson-Moskowitz if period not available from readers
        # WMO guide to wave analysis and forecasting pp. 14, WMO (1998)
        windspeed = self.wind_speed()
        # fallback value if no wind speed or Hs to avoid division by zero
        omega = 5*np.ones(windspeed.shape)
        omega[windspeed>0] = 0.877*9.81/(1.17*windspeed[windspeed>0])
        return omega

    def wave_period(self):
        if hasattr(self.environment, 'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment'
                ) and self.environment.sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment.max() > 0:
            # prefer using Tm02:
            T = self.environment.sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment.copy()
            logger.debug('Using mean period Tm02 as wave period')
        elif hasattr(self.environment, 'sea_surface_wave_period_at_variance_spectral_density_maximum'
                ) and self.environment.sea_surface_wave_period_at_variance_spectral_density_maximum.max() > 0:
            # alternatively use Tp
            T = self.environment.sea_surface_wave_period_at_variance_spectral_density_maximum.copy()
            logger.debug('Using peak period Tp as wave period')
        else:
            # calculate Tp from wind speed:
            logger.debug('Calculating wave period Tm02 from wind')
            T = (2*np.pi)/self._wave_frequency()
            self.environment.sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment = T

        #print '\n T %s \n' % str(T.mean())
        if T.min() == 0:
            logger.warning('Zero wave period found - '
                            'replacing with mean')
            T[T==0] = np.mean(T[T>0])

        logger.debug('   min: %f, mean: %f, max: %f' % (T.min(), T.mean(), T.max()))

        return T

    def wave_energy(self):
        return 9.81*1028*np.power(self.significant_wave_height(), 2)/16

    def wave_energy_dissipation(self):
        # Delvigne and Sweeney
        return 0.0034*self.sea_water_density()*9.81 * \
            np.power(self.significant_wave_height(), 2)

    def wave_damping_coefficient(self):
        omega = 2*np.pi / self.wave_period()
        return (10E-5)*omega * \
            np.power(self.wave_energy(), 0.25)

    # def sea_water_density(self):
    #    return 1027  # kg/m3

    def sea_surface_wave_breaking_fraction(self):
        # TODO: We should also have an option here for
        # the case when wave height is given, but no wind
        f = 0.032*(self.wind_speed() - 5)/self.wave_period()
        f[f < 0] = 0
        return f

    def air_density(self):
        return 1.225  # Could calculate from temperature

    def windspeed_from_stress(self):
        wind_stress = np.sqrt(
            self.environment.surface_downward_x_stress**2 +
            self.environment.surface_downward_y_stress**2)
        return windspeed_from_stress_polyfit(wind_stress)

    def solar_elevation(self):
        '''Solar elevation at present time and position of active elements.'''
        return solar_elevation(self.time, self.elements.lon, self.elements.lat)

    def sea_floor_depth(self):
        '''Sea floor depth (positive) for presently active elements'''

        if hasattr(self, 'environment') and \
                hasattr(self.environment, 'sea_floor_depth_below_sea_level'):
            if len(self.environment.sea_floor_depth_below_sea_level) == \
                    self.num_elements_active():
                sea_floor_depth = \
                    self.environment.sea_floor_depth_below_sea_level
        if 'sea_floor_depth' not in locals():
            env, env_profiles, missing = \
                self.env.get_environment(['sea_floor_depth_below_sea_level'],
                                     time=self.time, lon=self.elements.lon,
                                     lat=self.elements.lat,
                                     z=0*self.elements.lon, profiles=None)
            sea_floor_depth = \
                env['sea_floor_depth_below_sea_level'].astype('float32')
        return sea_floor_depth

    def sea_surface_height(self):
        '''Sea surface height (positive/negative) for presently active elements'''

        if hasattr(self, 'environment') and \
                hasattr(self.environment, 'sea_surface_height'):
            if len(self.environment.sea_surface_height) == \
                    self.num_elements_active():
                sea_surface_height = \
                    self.environment.sea_surface_height
        if 'sea_surface_height' not in locals():
            env, env_profiles, missing = \
                self.env.get_environment(['sea_surface_height'],
                                     time=self.time, lon=self.elements.lon,
                                     lat=self.elements.lat,
                                     z=0*self.elements.lon, profiles=None)
            sea_surface_height = \
                env['sea_surface_height'].astype('float32')
        return sea_surface_height


def wind_drag_coefficient(windspeed):
    '''Large and Pond (1981), J. Phys. Oceanog., 11, 324-336.'''
    Cd = 0.0012*np.ones(len(windspeed))
    Cd[windspeed > 11] = 0.001*(0.49 + 0.065*windspeed[windspeed > 11])
    return Cd

def windspeed_from_stress_polyfit(wind_stress):
    '''Inverting Large and Pond (1981) using polyfit'''
    windspeed = np.linspace(0, 30, 30)
    rho_air = 1.225
    stress = wind_drag_coefficient(windspeed)*rho_air*(windspeed**2)
    z = np.polyfit(stress, windspeed, 3)
    p = np.poly1d(z)
    return p(wind_stress)


def declination(time):
    '''Solar declination in degrees.'''
    try:
        day_of_year = time.timetuple().tm_yday
    except:
        day_of_year = np.asarray([t.timetuple().tm_yday for t in time])
    declination = \
        np.arcsin(np.deg2rad(-23.44)*
                  np.cos(np.radians((360.0/365.24)*(day_of_year + 10) +
                                    (360.0/np.pi)*0.0167*
                                    np.sin(np.radians((360.0/365.24)*
                                           (day_of_year - 2)))
                                    )))
    return np.rad2deg(declination)


def equation_of_time(time):
    '''Equation of time in minutes.'''
    time = np.atleast_1d(time)
    day_of_year = np.asarray([t.timetuple().tm_yday for t in time])
    hour = np.asarray([t.hour for t in time])
    gamma = 2*np.pi/365.0*(day_of_year - 1. +
                           (hour - 12.) / 24.)
    eqtime = 229.18 * (0.000075 + 0.001868*np.cos(gamma) -
                       0.032077*np.sin(gamma) -
                       0.014615*np.cos(2*gamma) - 0.040849*np.sin(2*gamma))
    return eqtime


def hour_angle(time, longitude):
    '''Solar hour angle in degrees.'''
    time = np.atleast_1d(time)
    time_offset = equation_of_time(time) + 4*longitude
    day_minutes = [t.hour*60.0 + t.minute + t.second/60.0 for t in time]
    true_solar_time = day_minutes + time_offset
    hour_angle = (true_solar_time/4.0) - 180.0  # degrees
    return hour_angle


def solar_elevation(time, longitude, latitude):
    '''Solar elevation in degrees.'''
    d_rad = np.deg2rad(declination(time))
    h = hour_angle(time, longitude)
    solar_elevation = np.rad2deg(np.arcsin(
        np.sin(np.deg2rad(latitude))*np.sin(d_rad) +
        np.cos(np.deg2rad(latitude))*np.cos(d_rad)*np.cos(np.deg2rad(h))))
    return solar_elevation
