import logging; logger = logging.getLogger(__name__)
import numpy as np
import cartopy.crs as ccrs
import pyproj

from . import srs

class RegularGrid:
  r"""
  A regular grid.

  The grid is projected using UTM with units in meters.

  TODO: Subclass ndarray so that we don't need to grid.grid.

  Attributes:

    grid: 2D grid, first dimension is easting (longitude), second northing (latitude).
  """
  srs  = None # pyproj.Proj
  crs  = None # pyproj.CRS
  ccrs = None # cartopy.crs.CRS
  grid = None               # values at grid, 2D np.array
  lon0, lat0 = None, None   # lower left corner
  x, y       = None, None   # x and y coordinates (1D)
  extent_xy  = None         # x and y corners.
  res = None                # resolution in meters

  # 2D arrays of longitudes and latitudes at gridcells.
  lons, lats = None, None

  def __init__(self, srs, grid):
    self.srs = srs
    self.crs = srs.crs
    self.grid = grid

  @staticmethod
  def new(lon, lat, res, shape):
    """
    Make a new regular grid.

    Args:

      lon, lat: Upper left corner of grid.

      res: Size of grid cell (meters).
    """
    assert len(shape) == 2

    p = srs.find_utm_proj(lon, lat)
    logger.info('Grid SRS: %s, resolution: %.2f m, upper left corner: lon=%.2f, lat=%.2f' % (p, res, lon, lat))

    r = RegularGrid(p, np.zeros(shape))
    r.lon0 = lon
    r.lat0 = lat
    r.res = res
    r.ccrs = srs.find_utm_ccrs(lon, lat)

    r.__make_grid__()

    return r

  def __make_grid__(self):
    x0, y0 = self.srs(self.lon0, self.lat0)
    x1 = x0 + self.grid.shape[0] * self.res
    y1 = y0 + self.grid.shape[1] * self.res
    self.extent_xy = [x0, x1, y0, y1]

    logger.debug("Grid [meters]: X = [%.1f ... %.1f], Y = [%.1f ... %.1f], shape = %s" % (x0, x1, y0, y1, self.grid.shape))

    self.x = np.linspace(x0, x1, self.grid.shape[0])
    self.y = np.linspace(y0, y1, self.grid.shape[1])

    xx, yy = np.meshgrid(self.x, self.y)
    self.lons, self.lats = self.srs(xx, yy, inverse = True)

  def center(self):
    """
    Center of grid in longitude and latitude
    """

    # find center of grid
    xc, yc = self.grid.shape
    xc, yc = self.x[xc//2], self.y[yc//2]

    return self.srs(xc, yc, inverse = True)

  def contains(self, x, y):
    """
    Check if points x and y are within grid.
    """
    x0, x1, y0, y1 = self.extent_xy

    return np.all(np.logical_and.reduce((x >= x0, x <= x1, y >= y0, y <= y1)))

  def plot(self, ax = None, crs = None):
    """
    Show grid on axis (or create new figure and axis)

    Args:

      ax: Axis to plot grid on.

      crs: Projection to use for plot (when setting up new plot), default (and
           fastest) is grid projection.
    """
    import matplotlib.pyplot as plt
    plt.style.use('dark_background')

    if ax is None:
      if crs is None:
        crs = self.ccrs

    ax = plt.axes(projection = crs, label = "%s" % np.random.randint(1000))
    im = ax.imshow(self.grid, origin = 'lower', transform = self.ccrs, extent = self.extent_xy, cmap = 'inferno')
    plt.colorbar(im, ax=ax, orientation = 'horizontal')
    ax.coastlines()
    ax.gridlines(draw_labels = True)

