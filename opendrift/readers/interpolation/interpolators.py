import logging
import numpy as np
from scipy.ndimage import map_coordinates, grey_dilation
import logging; logging.captureWarnings(True); logger = logging.getLogger(__name__)
from scipy.interpolate import interp1d, LinearNDInterpolator

logger = logging.getLogger('opendrift')  # using common logger

def expand_numpy_array(data):
    if isinstance(data, np.ma.MaskedArray):
        logger.warning('Converting masked array to numpy array before interpolating')
        data = np.ma.filled(data, fill_value=np.nan)
    if not np.isfinite(data).any():
        logger.warning('Only NaNs, returning')
        return
    mask = ~np.isfinite(data)
    data[mask] = np.finfo(np.float64).min
    data[mask] = grey_dilation(data, size=3)[mask]
    data[data==np.finfo(np.float64).min] = np.nan


###########################
# 2D interpolator classes
###########################

class Nearest2DInterpolator():

    def __init__(self, xgrid, ygrid, x, y):
        self.x = x
        self.y = y
        self.xi = (x - xgrid.min())/(xgrid.max()-xgrid.min())*len(xgrid)
        self.yi = (y - ygrid.min())/(ygrid.max()-ygrid.min())*len(ygrid)
        self.xi = np.round(self.xi).astype(np.uint32)
        self.yi = np.round(self.yi).astype(np.uint32)
        self.xi[self.xi >= len(xgrid)] = len(xgrid)-1
        self.yi[self.yi >= len(ygrid)] = len(ygrid)-1

    def __call__(self, array2d):
        return array2d[self.yi, self.xi]


class NDImage2DInterpolator():

    def __init__(self, xgrid, ygrid, x, y):
        self.x = x
        self.y = y
        self.xi = (x - xgrid.min())/(xgrid.max()-xgrid.min())*len(xgrid)
        self.yi = (y - ygrid.min())/(ygrid.max()-ygrid.min())*len(ygrid)

    def __call__(self, array2d):
        try:
            array2d = np.ma.array(array2d, mask=array2d.mask)
            array2d[array2d.mask] = np.nan  # Gives holes
        except:
            pass
        return np.ma.masked_invalid(
            map_coordinates(array2d, [self.yi, self.xi],
                            cval=np.nan, order=0))


class LinearND2DInterpolator():


    logger = logging.getLogger('opendrift')

    def __init__(self, xgrid, ygrid, x, y):
        self.block_x, self.block_y = np.meshgrid(xgrid, ygrid)
        self.block_x = self.block_x.ravel()
        self.block_y = self.block_y.ravel()
        self.x = x
        self.y = y

    def __call__(self, array2d):
        array_ravel = array2d.ravel()
        valid = np.isfinite(array_ravel)
        #if isinstance(array2d.mask, np.ndarray):
        #    valid = ~array2d.ravel().mask
        #elif array2d.mask == False:
        #    valid = np.ones(array_ravel.shape, dtype=bool)
        #elif array2d.mask == True:
        #    valid = np.zeros(array_ravel.shape, dtype=bool)

        if hasattr(self, 'interpolator'):
            if not np.array_equal(valid, self.interpolator.valid):
                logger.debug('Cannot reuse interpolator - validity of '
                              'array is different from original.')
        if hasattr(self, 'interpolator') and (np.array_equal(
                 valid, self.interpolator.valid)):
            # Reuse stored interpolator with new data
            self.interpolator.values[:, 0] = \
                (array_ravel[valid])
        else:
            # Make new interpolator for given x,y
            self.interpolator = LinearNDInterpolator(
                (self.block_y[valid],
                 self.block_x[valid]),
                array_ravel[valid])
            # Store valid array, to determine if can be used again
            self.interpolator.valid = valid
            # Call interpolator to avoid threading-problem:
            # https://github.com/scipy/scipy/issues/8856
            self.interpolator((0,0))

        return self.interpolator(self.y, self.x)


class Linear2DInterpolator():

    logger = logging.getLogger('opendrift')

    def __init__(self, xgrid, ygrid, x, y):
        self.x = x
        self.y = y
        self.xi = (x - xgrid[0])/(xgrid[-1]-xgrid[0])*(len(xgrid)-1)
        self.yi = (y - ygrid[0])/(ygrid[-1]-ygrid[0])*(len(ygrid)-1)

    def __call__(self, array2d):
        if isinstance(array2d,np.ma.MaskedArray):
            logger.debug('Converting masked array to numpy array for interpolation')
            array2d = np.ma.filled(array2d, fill_value=np.nan)
        if not np.isfinite(array2d).any():
            logger.warning('Only NaNs input to linearNDFast - returning')
            return np.nan*np.ones(len(self.xi))

        # Fill NaN-values with nearby real values
        interp = map_coordinates(array2d, [self.yi, self.xi],
                                 cval=np.nan, order=1)
        missing = np.where(~np.isfinite(interp))[0]
        i=0
        while len(missing) > 0:
            i += 1
            if i > 10:
                logger.warning('Still NaN-values after 10 iterations, exiting!')
                return interp
            logger.debug('Linear2DInterpolator informational: NaN values for %i elements, expanding data %i' %
                          (len(missing), i))
            expand_numpy_array(array2d)
            interp[missing] = map_coordinates(
                array2d, [self.yi[missing], self.xi[missing]],
                cval=np.nan, order=1, mode='nearest')
            missing = np.where(~np.isfinite(interp))[0]

        return interp


horizontal_interpolation_methods = {
    'nearest': Nearest2DInterpolator,
    'ndimage': NDImage2DInterpolator,
    'linearND': LinearND2DInterpolator,
    'linearNDFast': Linear2DInterpolator}


###########################
# 1D interpolator classes
###########################

class Nearest1DInterpolator():

    def __init__(self, zgrid, z):

        # Truncating above and below
        z[z < zgrid.min()] = zgrid.min()
        z[z > zgrid.max()] = zgrid.max()
        # Interpolator zgrid -> index
        if zgrid[1] > zgrid[0]:  # increasing
            z_interpolator = interp1d(zgrid, range(len(zgrid)))
        else:  # decreasing values, must flip for interpolator
            z_interpolator = interp1d(zgrid[::-1], range(len(zgrid))[::-1])
        # Indices corresponding to nearest value in zgrid
        self.zi = np.round(z_interpolator(z)).astype(np.uint8)
        self.zi[self.zi < 0] = 0
        self.zi[self.zi >= len(zgrid)] = len(zgrid) - 1

    def __call__(self, array2d):
        return array2d[self.zi, range(len(self.zi))]


class Linear1DInterpolator():

    def __init__(self, zgrid, z):

        # Truncating above and below
        z[z < zgrid.min()] = zgrid.min()
        z[z > zgrid.max()] = zgrid.max()
        # Interpolator zgrid -> index
        if zgrid[1] > zgrid[0]:  # increasing
            z_interpolator = interp1d(zgrid, range(len(zgrid)))
        else:  # decreasing values, must flip for interpolator
            z_interpolator = interp1d(zgrid[::-1], range(len(zgrid))[::-1])
        z_interpolator(z[0])  # to prevent threading issues
        # Indices corresponding to layers above and below
        interp_zi = z_interpolator(z)
        self.index_above = np.floor(interp_zi).astype(np.int8)
        self.index_above[self.index_above < 0] = 0
        self.index_below = np.minimum(self.index_above + 1, len(zgrid) - 1)
        self.weight_above = 1 - (interp_zi - self.index_above)
        self.xi = range(len(z))

    def __call__(self, array2d):
        return array2d[self.index_above, self.xi]*self.weight_above + \
               array2d[self.index_below, self.xi]*(1 - self.weight_above)

vertical_interpolation_methods = {
    'nearest': Nearest1DInterpolator,
    'linear': Linear1DInterpolator}


def fill_NaN_towards_seafloor(array):
    """Extrapolate NaN-values (missing) towards seafloor"""
    filled = False
    for i in range(1, array.shape[0]):
        mask = np.isnan(array[i,:,:])
        if np.sum(mask) > 0:
            array[i, mask] = array[i-1, mask]
            filled = True
    return filled

