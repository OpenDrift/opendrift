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

from __future__ import (absolute_import, division)

import numpy as np


def sdepth(depth, sigma, stagger="rho", Vtransform=1):
    """Depth of s-levels

    *H* : arraylike
      Bottom depths [meter, positive]

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

    >>> fid = Dataset(roms_file)
    >>> H = fid.variables['h'][:, :]
    >>> C = fid.variables['Cs_r'][:]
    >>> Hc = fid.variables['hc'].getValue()
    >>> z_rho = sdepth(H, Hc, C)

    """
    H = np.asarray(depth)
    H_shape = H.shape      # Save the shape of H
    H = depth.ravel()         # and make H 1D for easy shape maniplation
    C = np.asarray(sigma)
    N = len(C)
    outshape = (N,) + H_shape       # Shape of output
    if stagger == 'rho':
        S = -1.0 + (0.5+np.arange(N))/N    # Unstretched coordinates
    elif stagger == 'w':
        S = np.linspace(-1.0, 0.0, N)
    else:
        raise ValueError("stagger must be 'rho' or 'w'")

    if Vtransform == 1:         # Default transform by Song and Haidvogel
        A = (S - C)[:, None]
        B = np.outer(C, H)
        return (A + B).reshape(outshape)

    elif Vtransform == 2:       # New transform by Shchepetkin
        N = S[:, None] + np.outer(C, H)
        D = (1.0 + 1/H)
        return (N/D).reshape(outshape)

    else:
        raise ValueError("Unknown Vtransform")


# ------------------------------------------
# Vertical slicing e.t.c.
# ------------------------------------------

# -----------------------------------------------


def multi_zslice(F, S, Z):

    """Slice a 3D ECOM field to fixed depth

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

