from abc import abstractmethod
from numbers import Number
from typing import List
import pyproj
import xarray as xr
import numpy as np


class Combine:
    """Combining two readers into a third one. You can use usual operators,
    but also more complex ones such as gaussian combining.
    """

    def __add__(self, other):
        from .readerops import Combined as ReaderCombined
        from .numops import Combined as NumCombined
        from ..basereader import BaseReader

        if isinstance(other, Number):
            return NumCombined.add(other, self)
        elif isinstance(other, BaseReader):
            return ReaderCombined(self, other, lambda a, b: a + b)
        else:
            return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        from .numops import Combined as NumCombined
        if isinstance(other, Number):
            return NumCombined.mul(other, self)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        from .numops import Combined as NumCombined

        if isinstance(other, Number):
            return NumCombined.div(other, self)
        else:
            return NotImplemented

    def __sub__(self, other):
        return self + (-1 * other)

    def combine_gaussian(self, measurement_reader, std):
        """Mix two readers with a gaussian, whose std is the one given as an argument.
        The measurment reader have to be of type timeseries, with a lon and lat
        attributes that are taken as the center of the measure.
        """
        from .readerops import Combined as ReaderCombined

        def gaussian_factor(lon, lat, lon_center, lat_center, std):
            geod = pyproj.Geod(ellps='WGS84')

            assert isinstance(np.broadcast_arrays(lon, lat), list), f"requested lon and lat not broadcastable"
            lon, lat = np.broadcast_arrays(lon, lat)
            requested_shape = lon.shape
            requested_ndim = len(requested_shape)

            ##
            if isinstance(lon_center, float) :
                    lon_center = lon_center * np.ones(requested_shape)
            #
            elif lon_center.shape != requested_shape :
                lon_center = np.expand_dims(lon_center,  tuple(range(1, requested_ndim+1)))
                lon, lon_center = np.broadcast_arrays(lon_center, lon)

            ##
            if isinstance(lat_center, float) :
                    lat_center = lat_center * np.ones(requested_shape)
            #
            elif lat_center.shape != requested_shape :
                lat_center = np.expand_dims(lat_center,  tuple(range(1, requested_ndim+1)))
                lat, lat_center = np.broadcast_arrays(lat_center, lat)

            ##
            if isinstance(std, float) :
                std = std * np.ones(requested_shape)
            elif lon_center.shape != requested_shape :
                std = np.expand_dims(std,  tuple(range(1, requested_ndim+1)))
                std, _ = np.broadcast_arrays(std, lat)

            _, _, distance = geod.inv(lon, lat, lon_center, lat_center)
            exponential_factor = np.exp( -np.power(distance/std, 2.) / 2)
            return exponential_factor

        return ReaderCombined(a = self, b = measurement_reader, op = gaussian_factor, op_type= "combine_gaussian", external_params = std)



class Filter:
    @property
    @abstractmethod
    def variables(self) -> List[str]:
        pass

    def filter_vars(self, vars):
        """
        Only keep the specified variables.
        """
        from .filter import FilterVariables
        return FilterVariables(self, vars)

    def exclude_vars(self, vars):
        """
        Remove the specified variables.
        """
        from .filter import FilterVariables

        vars = list(set(self.variables) - set(vars))
        return FilterVariables(self, vars)

