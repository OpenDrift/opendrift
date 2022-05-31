import numpy as np


def vec_nearest(X, xp):
    """
    Find nearest element in vector `X` to `xp`.

    Args:

        X   Array
        xp  Array

    Returns:
        i   M   vector of indices [0..N] of closest element in
                X[0..N, i] to xp[i]
    """
    xp = np.atleast_1d(xp)
    return np.argmin(np.abs(X - xp), axis=0)

