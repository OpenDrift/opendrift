import logging
logger = logging.getLogger(__name__)
from .grid import RegularGrid
from . import srs
from .interp import vec_nearest
from abc import abstractmethod
import numpy as np
from datetime import timedelta, datetime


class Simulation:
    r"""
    The simulation.

    It contains the problem to be simulated, means to read necessary input
    variables, and the physics for modeling the convection of the initial
    conditions.

    Convection:

      Convection consists of advection and diffusion.

      Diffusion is given by:

      .. math::

        \frac{\partial U}{\partial t} = D \left( \frac{\partial^2 U}{\partial x^2} + \frac{\partial^2 U}{\partial y^2} \right)


      The convection equation is (`wiki
      <https://en.wikipedia.org/wiki/Convection%E2%80%93diffusion_equation>`_):

      .. math::

        \frac{\partial c}{\partial t} = ...

      with the assumptions that:

        * the diffusion constant `D` is constant for the field,
        * and that the flow `u` is incompressible (i.e. has *no divergence*).

      the equation simplifies to:

      .. math::

        \frac{\partial c}{\partial t} = D \nabla^2 c - \mathbf{v} \cdot \nabla T

      where :math:`\nabla^2 = \triangle` is the Laplacian.
    """

    """The model grid."""
    grid = None

    """Reference time `datetime.Datetime`"""
    t0 = None

    """Current time offset after `t0`"""
    t = 0.

    """
    Diffusivity (:math:`m^2/s`). E.g. between 0.01 and 0.1 for oil on the
    surface of the ocean (`Matsuzakia et. al., 2017
    <https://www.sciencedirect.com/science/article/pii/S0025326X16308426>`_).

    Decreasing diffusivity places stricter stability criteria on time step.
    """
    D = 0.1

    """Porosity, rate of liquid volume to total volume (fraction of flux)"""
    rho = 1.

    """List of readers"""
    readers = None

    def __init__(self, grid):
        self.grid = grid
        self.readers = []
        self.t0 = datetime.utcnow()

    @classmethod
    def new(cls, lon0 = 10., lat0 = 65., res = 10., shape=(100, 100)):
        """
        New simulation on a :class:`.grid.RegularGrid`.

        Args:

            lon0: lower-left corner longitude

            lat0: lower-left corner latitude

            res: resolution (dx and dy)

            shape: shape (size) of grid
        """
        s = cls(grid=RegularGrid.new(lon0, lat0, res, shape))
        logger.info("New %s with Regular Grid" % (cls.__name__))

        return s

    def source(self, lon, lat, X):
        """
        Source `X` onto grid with lower-left corner `lon`, `lat`.
        """
        x0, y0 = self.grid.srs(lon, lat)
        assert self.grid.contains(x0, y0)

        ix0 = vec_nearest(self.grid.x, x0)
        iy0 = vec_nearest(self.grid.y, y0)
        print(ix0)
        ix1 = ix0 + X.shape[0]
        iy1 = iy0 + X.shape[1]

        assert ix1 < self.grid.grid.shape[0]
        assert iy1 < self.grid.grid.shape[1]

        self.grid.grid[ix0:ix1, iy0:iy1] = X

    def source_gaussian_blob(self, lon, lat, A=1., N=10, sigma=10.):
        r"""
        Source a Gaussian blob (2D normal distribution) at `lon` and `lat` with
        `sigma` (standard deviation, meters) radius.

        Args:

          lon, lat: Center coordinates, or :math:`\bar \mu`.

          A: Amplitude.

          N: Kernel size.

          sigma: standard deviation (:math:`\sigma`) in meters.
        """
        def gaussian2d(n, std):
            """
            Gaussian 2D window with size `n` and standard deviation `std`.
            """

            from scipy import signal
            gk = signal.windows.gaussian(n, std=std)
            return np.outer(gk, gk)

        # center of kernel
        x, y = self.grid.srs(lon, lat)
        assert self.grid.contains(x, y)

        ix0 = vec_nearest(self.grid.x, x) - N // 2
        iy0 = vec_nearest(self.grid.y, y) - N // 2
        ix1 = ix0 + N
        iy1 = iy0 + N

        assert ix0 > 0 and iy0 > 0

        S = gaussian2d(N, sigma / self.grid.res)

        self.grid.grid[ix0:ix1, iy0:iy1] = S

    def U(self, t):
        """
        Get U (ocean current) for t0 + t
        """
        Ux, Uy = None, None
        t = self.t0 + timedelta(seconds = t)

        for r in self.readers:
            vv = r.variables()
            if 'x_sea_water_velocity' in vv and 'y_sea_water_velocity' in vv:
                Ux, Uy = r.read_grid(self.grid, ['x_sea_water_velocity', 'y_sea_water_velocity'], t)

            if Ux is not None and Uy is not None:
                break

        return Ux, Uy

    @abstractmethod
    def step(self, dt=None):
        """
        Step the simulation.

        Stepping the simulation involves applying diffusion and advection to the field.

        Args:

          dt: time delta (or use automatic).
        """
        pass

    @abstractmethod
    def integrate(self, dt=None, max_t=None, max_steps=None, observer = None):
        """
        Run simulation until termination condition is met.

        Args:

          dt: override time step

          observer: function to call after each step taking simulation object as first argument. the function may return False to stop the integration.

          Termination conditions:

          max_t: max time
          max_steps: max iterations
        """

        if max_t is not None and max_steps is not None:
            logger.warning(
                "no termination condition supplied, using max_steps = 1000")
            max_steps = 1000

        logger.info("integrating.. (dt = %s)", dt)

        n = 0

        while True:
            logger.debug("step: n = %d, t = %s (dt = %s)" % (n, self.t, dt))

            if max_t is not None and self.t >= max_t:
                logger.info("max_t reached.")
                break

            if max_steps is not None and n >= max_steps:
                logger.info("max_steps reached.")
                break

            self.step(dt=dt)

            n += 1


class ExplSimulation(Simulation):
    r"""
    A simple explicit scheme for integrating the convection-equation.

    * Forward difference in time
    * `ndimage.laplace` and `np.gradient` for spatial differences.

    https://en.wikipedia.org/wiki/Numerical_solution_of_the_convection%E2%80%93diffusion_equation#Solving_the_convection%E2%80%93diffusion_equation_using_the_finite_difference_method

    .. seealso::

      :class:`Simulation`.

    """
    def stability(self, dx, D, umax):
        """
        https://en.wikipedia.org/wiki/Numerical_solution_of_the_convection%E2%80%93diffusion_equation#Solving_the_convection%E2%80%93diffusion_equation_using_the_finite_difference_method
        """
        h = 2 * D / (self.rho * umax)
        dt = dx**2 / (2 * D)

        return h, dt

    def step(self, dt=None):
        from scipy import ndimage

        dx, dy = self.grid.res, self.grid.res
        D = self.D
        # Ux = np.sqrt(50 * 1.) * np.ones(self.grid.grid.shape)
        # Uy = Ux
        Ux, Uy = self.U(self.t)
        maxu = np.max(np.sqrt(Ux**2 + Uy**2).ravel())
        logger.debug("maxu = %s" % maxu)

        h, ddt = self.stability(dx, D, maxu)

        if h < dx:
            logger.warning("dx too big, dx = %s > h = %s" % (dx, h))

        if dt is None:
            dt = ddt
        elif ddt < dt:
            logger.warning("dt too big, dt = %s > ddt = %s" % (dt, ddt))

        # diffusion
        diff = D * 1. / dx**2 * ndimage.laplace(self.grid.grid)

        # advection
        gfx, gfy = np.gradient(self.grid.grid, dx, dy)
        adv = -(gfx * Ux + gfy * Uy)

        self.grid.grid = self.grid.grid + dt * (diff + adv)

        self.t += dt
