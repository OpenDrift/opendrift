from .state import State
import xarray as xr


class Result(State):
    data: xr.Dataset

    def __init__(self, data):
        self.data = data

    @staticmethod
    def from_simulation(path):
        pass

    @staticmethod
    def from_dataset(ds: xr.Dataset):
        """

        You can open a file with:

        >>> ds = xr.open_dataset(f)
        >>> r = Result.from_dataset(ds)
        """
        return Result(ds)

    @staticmethod
    def from_file(path):
        """
        Open a OpenDrift dataset from file.

        >>> f = 'opendrift.nc'
        >>> r = Result.from_file(f)
        """
        return Result.from_dataset(xr.open_dataset(path))

    def to_file(self, path):
        """
        Save result to file.
        """
        return self.data.to_file(path, format='netCDF4')

