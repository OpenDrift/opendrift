from abc import abstractmethod
import numpy as np
import logging
logger = logging.getLogger(__name__)


class Reader:
    @abstractmethod
    def variables(self):
        """
        Get list of variables
        """
        return []

    @abstractmethod
    def read_grid(self, grid, var, t):
        """
        Read variable for grid.

        Args:

            grid: grid to read variable for

            var: list of strs, name of variables to read

            t: datetime (UTC)

        Returns:

            v0, v1, ...

            tuple with each variables specified in var with same shape as `grid`.
        """
        pass


class ConstantReader(Reader):
    def __init__(self, consts):
        """
        Args:

            consts: dict, name of var: constant
        """
        logger.info("constant reader for: %s" % ", ".join(consts.keys()))
        self.consts = consts
        super().__init__()

    def read_grid(self, grid, var, _):
        for v in var:
            assert v in self.consts, "missing variable %s" % v

        return tuple(self.consts[v] * np.ones(grid.grid.shape) for v in var)

    def variables(self):
        return list(self.consts.keys())

    @staticmethod
    def new_xy(x=.5, y=.5):
        return ConstantReader({
            'x_sea_water_velocity': x,
            'y_sea_water_velocity': y
        })


class OpendriftReader(Reader):
    """
    Wrapper around an Opendrift reader.
    """
    def __init__(self, reader):
        """
        Args:

            reader: opendrift reader
        """
        logger.info("opendrift reader: %s (%s)" %
                    (reader.name, ", ".join(reader.variables)))

        self.r = reader

        super().__init__()

    def variables(self):
        return self.r.variables

    def read_grid(self, grid, var, t):
        logger.debug('reading %s for %s' % (var, t))

        x, y = self.r.lonlat2xy(grid.lons.ravel(), grid.lats.ravel())

        env, _ = self.r.get_variables_interpolated_xy(
            variables=var,
            time=t,
            x=x,
            y=y,
            z=np.zeros(grid.grid.shape).ravel(),
            rotate_to_proj=grid.srs)

        u = tuple(env[v].reshape(grid.grid.shape).filled(fill_value=np.nan) for v in var)

        for uu, vv in zip(u, var):
            if np.isnan(uu).any():
                logger.warning("nan's in %s" % vv)

        return u
