"""
Readers may be combined by using operators and the helper functions.
"""

from ..basereader import BaseReader

def mean(*readers, weights = None, variables = None) -> BaseReader:
    """
    Combine readers and calculate the mean from the listed variables (or all if not specified).

    Arguments:
        readers: An iterables of readers, _taken by reference_.

        weights: A list of weights to apply to readers. Will be normalized.

        variables: Which variables to include in this mean, e.g.: `[x_wind,
        y_wind]`, if not given the intersection of available variables is used.

    Returns: A reader implementing :class:`BaseReader`.
    """
    raise NotImplemented

