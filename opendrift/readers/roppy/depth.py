# -*- coding: utf-8 -*-

"""Vertical structure functions for ROMS

:func:`sdepth`
  Depth of s-levels
:func:`zslice`
  Slice a 3D field in s-coordinates to fixed depth
:func:`multi_zslice`
  Slice a 3D field to several depth levels
:func:`z_average`
  Vertical average of a 3D field
:func:`s_stretch`
  Compute vertical stretching arrays Cs_r or Cs_w

"""

# -----------------------------------
# Bjørn Ådlandsvik <bjorn@imr.no>
# Institute of Marine Research
# Bergen, Norway
# 2010-09-30
# -----------------------------------

import numpy as np

def sdepth(H, zeta, Hc, C, stagger="rho", Vtransform=1):
    """Depth of s-levels

    *H* : arraylike
      Bottom depths [meter, positive]
    
    *zeta* : scalar, arraylike
      Surface elevation [meter]

    *Hc* : scalar
       Critical depth

    *cs_r* : 1D array
       s-level stretching curve

    *stagger* : [ 'rho' | 'w' ]

    *Vtransform* : [ 1 | 2 ]
       defines the transform used, defaults 1 = Song-Haidvogel

    Returns an array with ndim = H.ndim + 1 and
    shape = cs_r.shape + H.shape with the depths of the
    mid-points in the s-levels.

    Typical usage::

    .. code::

        fid = Dataset(roms_file)
        H = fid.variables['h'][:, :]
        zeta = fid.variables['zeta'][:, :]
        C = fid.variables['Cs_r'][:]
        Hc = fid.variables['hc'].getValue()
        z_rho = sdepth(H, zeta, Hc, C)

    """
    H = np.asarray(H)
    zeta = np.asarray(zeta)
    Hshape = H.shape      # Save the shape of H
    Hflat = H.ravel()         # and make H 1D for easy shape manipulation
    C = np.asarray(C)
    N = len(C)
    outshape = (N,) + Hshape       # Shape of output
    if stagger == 'rho':
        S = -1.0 + (0.5+np.arange(N))/N    # Unstretched coordinates
    elif stagger == 'w':
        S = np.linspace(-1.0, 0.0, N)
    else:
        raise ValueError("stagger must be 'rho' or 'w'")

    # https://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#_ocean_s_coordinate_generic_form_1
    if Vtransform == 1:         # Default transform by Song and Haidvogel
        A = Hc * (S - C)[:, None]
        B = np.outer(C, H)
        Zo_rho = (A + B).reshape(outshape)
        z_rho = Zo_rho + zeta * (1 + Zo_rho / H)
        return z_rho

    # https://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#_ocean_s_coordinate_generic_form_2
    elif Vtransform == 2:       # New transform by Shchepetkin
        N = Hc*S[:, None] + np.outer(C, Hflat)
        D = (Hflat + Hc)
        z_rho = zeta + (zeta + H) * (N/D).reshape(outshape)
        return z_rho

    else:
        raise ValueError("Unknown Vtransform")

# ------------------------------------


def sdepth_w(H, Hc, cs_w):
    """Return depth of w-points in s-levels

    Kept for backwards compatibility
    use *sdepth(H, Hc, cs_w, stagger='w')* instead

    """
    return sdepth(H, Hc, cs_w, stagger='w')

# ------------------------------------------
# Vertical slicing e.t.c.
# ------------------------------------------


def zslice(F, S, z):
    """Vertical slice of a 3D ROMS field

    Vertical interpolation of a field in s-coordinates to
    (possibly varying) depth level

    *F* : array with vertical profiles, first dimension is vertical

    *S* : array with depths of the F-values,

    *z* : Depth level(s) for output, scalar or ``shape = F.shape[1:]``
          The z values should be negative

    Return value : array, `shape = F.shape[1:]`, the vertical slice

    Example:
    H is an array of depths (positive values)
    Hc is the critical depth
    C is 1D containing the s-coordinate stretching at rho-points
    returns F50, interpolated values at 50 meter with F50.shape = H.shape

    .. code::

       z_rho = sdepth(H, Hc, C)
       F50 = zslice(F, z_rho, -50.0)

    """

    # TODO:
    # Option to Save A, D, Dm
    #   => faster interpolate more fields to same depth

    F = np.asarray(F)
    S = np.asarray(S)
    z = np.asarray(z, dtype='float')
    Fshape = F.shape  # Save original shape
    if S.shape != Fshape:
        raise ValueError("F and z_r must have same shape")
    if z.shape and z.shape != Fshape[1:]:
        raise ValueError("z must be scalar or have shape = F.shape[1:]")

    # Flatten all non-vertical dimensions
    N = F.shape[0]        # Length of vertical dimension
    M = F.size // N        # Combined length of horizontal dimension(s)
    F = F.reshape((N, M))
    S = S.reshape((N, M))
    if z.shape:
        z = z.reshape((M,))

    # Find integer array C with shape (M,)
    # with S[C[i]-1, i] < z <= S[C[i], i]
    # C = np.apply_along_axis(np.searchsorted, 0, S, z)
    # but the following is much faster
    C = np.sum(S < z, axis=0)
    C = C.clip(1, N-1)

    # For vectorisation
    # construct index array tuples D and Dm such that
    #   F[D][i]  = F[C[i], i]
    #   F[Dm][i] = F[C[i]-1, i]
    I = np.arange(M, dtype='int')
    D = (C, I)
    Dm = (C-1, I)

    # Compute interpolation weights
    A = (z - S[Dm]) / (S[D]-S[Dm])
    A = A.clip(0.0, 1.0)   # Control the extrapolation

    # Do the linear interpolation
    R = (1-A)*F[Dm]+A*F[D]

    # Give the result the correct s
    R = R.reshape(Fshape[1:])

    return R

# -----------------------------------------------


def multi_zslice(F, S, Z):

    """Slice a 3D ROMS field to fixed depth

    Vertical interpolation of a field in s-coordinates to
    fixed vertical level

    *F* : array of with vertical profiles, first dimension is vertical

    *S* : array with depth of s-levels (at rho-points)
        1D (constant depth) or  S.shape = F.shape

    *Z* : single depth value, negative

    Returns : array, ``shape = F.shape[1:]`` the vertical slice

    """

    # TODO:
    # Option to Save A, D, Dm
    #   => faster interpolate more fields to same depth

    F = np.asarray(F)
    S = np.asarray(S)
    Fshape = F.shape  # Save original shape

    # Flat all dimensions after first
    N = F.shape[0]
    M = F.size // N
    F = F.reshape((N, M))
    S = S.reshape((N, M))

    # Make z.shape = (M,)
    Z = np.asarray(Z, dtype='float')
    # Valid possibilities
    # 1) Z = single scalar (shape = ()), one constant value
    # 2) Z = 1D array, shape=(kmax), a set of constant depths
    # 3) Z = 2D or more, reshapeable to (kmax, M)

    if Z.ndim == 0:
        Z = Z + np.zeros((1, M))
        kmax = 1
    elif Z.ndim == 1:
        kmax = Z.size
        Z = Z[:, np.newaxis] + np.zeros((kmax, M))
    else:
        kmax = Z.size // M
        Z = Z.reshape((kmax, M))

    # Find C, C.shape = (kmax, M) such that
    # z_r[C[k,i]-1, i] < Z[k] <= z_r[C[k,i], i]

    # shape: kmax, N, M => kmax, M
    C = np.sum(S[np.newaxis, :, :] < Z[:, np.newaxis, :], axis=1)
    C = C.clip(1, N-1)

    # Horizontal index
    I = np.arange(M, dtype=int)

    # Compute interpolation weights
    A = (Z - S[(C-1, I)])/(S[(C, I)]-S[(C-1, I)])
    A = A.clip(0.0, 1.0)   # Control the extrapolation

    # Do the interpolation
    R = (1-A)*F[(C-1, I)]+A*F[(C, I)]

    # Give the result the correct shape
    R = R.reshape((kmax,) + Fshape[1:])

    return R, (A, C, I, kmax)

# ------------------------------------------------------


def z_average(F, z_r, z0, z1):
    """Slice a 3D ROMS field to fixed depth

    Vertical interpolation of a field in s-coordinates to
    fixed vertical level

    *F* : array
      Vertical profiles, first dimension is vertical
    *z_r* : array
      Depth of s-levels (at rho-points), requires `z_r.shape = F.shape`
    *z0*, *z1* : floats
      Single depth values with z0 <= z1 <= 0

    return value : array
      `shape = F.shape[1:]`, the vertical average

    """

    F = np.asarray(F)
    z_r = np.asarray(z_r)
    Fshape = F.shape  # Save original shape

    # Flatten all dimensions after first
    N = F.shape[0]
    M = F.size // N
    F = F.reshape((N, M))
    z_r = z_r.reshape((N, M))

    # z0, z1 are scalars or horizontal arrays
    z0 = np.asarray(z0)
    if z0.shape:  # Array, must be 2D
        z0 = z0.reshape((M,))
    z1 = np.asarray(z1)
    if z1.shape:
        z1 = z1.reshape((M,))

    # Bracket z0, i.e.
    # Find integer array C0 with shape (M,)
    # with z_r[C0[i]-1, i] < z0 <= z_r[C0[i], i]
    # Can be done with:
    #   C0 = np.apply_along_axis(np.searchsorted, 0, z_r, z0)
    # but the following is much faster
    C0 = np.sum(z_r < z0, axis=0)
    C0 = C0.clip(1, N-1)       # Clip to avoid illegal indices

    # Bracket z1
    C1 = np.sum(z_r < z1, axis=0)
    C1 = C1.clip(1, N-1)

    # Use advanced indexing for vectorisation
    #   F[(C0,I)][i]  = F[C0[i], i]
    I = np.arange(M, dtype='int')
    # Interpolate F to the two levels
    A0 = (z0 - z_r[(C0-1, I)]) / (z_r[(C0, I)]-z_r[(C0-1, I)])
    A0 = A0.clip(0.0, 1.0)   # Control the extrapolation
    F0 = (1-A0)*F[(C0-1, I)]+A0*F[(C0, I)]
    A1 = (z1 - z_r[(C1-1, I)])/(z_r[(C1, I)]-z_r[(C1-1, I)])
    A1 = A1.clip(0.0, 1.0)
    F1 = (1-A1)*F[(C1-1, I)] + A1*F[(C1, I)]

    # Find indices again (unclipped)
    C0 = np.sum(z_r < z0, axis=0)
    C1 = np.sum(z_r < z1, axis=0)

    R = np.zeros(M, dtype=np.float64)
    X = np.zeros(N+2, dtype=np.float64)
    Y = np.zeros(N+2, dtype=np.float64)
    z0 = z0 + R    # Make sure they are spatial arrays
    z1 = z1 + R    # For indexing below

    for i in I:
        X[:] = 0.0
        Y[:] = 0.0
        nz = C1[i] - C0[i]  # Number of rho-points between z0 and z1
        # Set up arrays for trapezoidal integration
        X[0] = z0[i]
        X[1:nz+1] = z_r[C0[i]:C1[i], i]
        X[nz+1] = z1[i]
        Y[0] = F0[i]
        Y[1:nz+1] = F[C0[i]:C1[i], i]
        Y[nz+1] = F1[i]
        # Perform the integration
        R[i] = 0.5 * np.dot(X[1:nz+2]-X[0:nz+1], Y[1:nz+2]+Y[0:nz+1])

    # Compute average and revert to correct shape
    R = R / (z1-z0)
    R = R.reshape(Fshape[1:])

    return R

# ----------------------------------


def s_stretch(N, theta_s, theta_b, stagger='rho', Vstretching=1):
    """Compute a s-level stretching array

    *N* : Number of vertical levels

    *theta_s* : Surface stretching factor

    *theta_b* : Bottom stretching factor

    *stagger* : "rho"|"w"

    *Vstretching* : 1|2|4

    """

    if stagger == 'rho':
        S = -1.0 + (0.5+np.arange(N))/N
    elif stagger == "w":
        S = np.linspace(-1.0, 0.0, N+1)
    else:
        raise ValueError("stagger must be 'rho' or 'w'")

    if Vstretching == 1:
        cff1 = 1.0 / np.sinh(theta_s)
        cff2 = 0.5 / np.tanh(0.5*theta_s)
        return ((1.0-theta_b)*cff1*np.sinh(theta_s*S)
                + theta_b*(cff2*np.tanh(theta_s*(S+0.5))-0.5))

    elif Vstretching == 2:
        a, b = 1.0, 1.0
        Csur = (1 - np.cosh(theta_s * S))/(np.cosh(theta_s) - 1)
        Cbot = np.sinh(theta_b * (S+1)) / np.sinh(theta_b) - 1
        mu = (S+1)**a * (1 + (a/b)*(1-(S+1)**b))
        return mu*Csur + (1-mu)*Cbot

    elif Vstretching == 4:
        C = (1 - np.cosh(theta_s * S)) / (np.cosh(theta_s) - 1)
        C = (np.exp(theta_b * C) - 1) / (1 - np.exp(-theta_b))
        return C

    else:
        raise ValueError("Unknown Vstretching")


# wrapper for backwards compatibility
def s_stretch_w(N, theta_s, theta_b, Vstretching=1):
    """Obsolete use *s_stretch* instead"""

    return s_stretch(N, theta_s, theta_b,
                     stagger='w', Vstretching=Vstretching)
