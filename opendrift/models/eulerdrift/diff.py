import numpy as np
from scipy import ndimage

# TODO: Check out `findiff` package.

def grad2d(f, dx, dy):
  return (
      ndimage.gaussian_filter1d(f, sigma = 1, axis = 0, order = 1, mode = 'constant', cval = 0.) / dx,
      ndimage.gaussian_filter1d(f, sigma = 1, axis = 1, order = 1, mode = 'constant', cval = 0.) / dy)

