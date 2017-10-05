import logging

import numpy as np
from scipy.ndimage import map_coordinates
import scipy.ndimage as ndimage
from scipy.interpolate import interp1d, LinearNDInterpolator


###########################
# 2D interpolator classes
###########################

class Nearest2DInterpolator():

    def __init__(self, xgrid, ygrid, x, y):
        self.x = x
        self.y = y
        self.xi = (x - xgrid.min())/(xgrid.max()-xgrid.min())*len(xgrid)
        self.yi = (y - ygrid.min())/(ygrid.max()-ygrid.min())*len(ygrid)
        self.xi = np.round(self.xi).astype(np.int)
        self.yi = np.round(self.yi).astype(np.int)
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

    def __init__(self, xgrid, ygrid, x, y):
        self.block_x, self.block_y = np.meshgrid(xgrid, ygrid)
        self.block_x = self.block_x.ravel()
        self.block_y = self.block_y.ravel()
        self.x = x
        self.y = y

    def __call__(self, array2d):
        array_ravel = array2d.ravel()
        if isinstance(array2d.mask, np.ndarray):
            valid = ~array2d.ravel().mask
        elif array2d.mask == False:
            valid = np.ones(array_ravel.shape, dtype=bool)
        elif array2d.mask == True:
            valid = np.zeros(array_ravel.shape, dtype=bool)

        if hasattr(self, 'interpolator'):
            if not np.array_equal(valid, self.interpolator.valid):
                logging.debug('Cannot reuse interpolator - validity of '
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

        return self.interpolator(self.y, self.x)
        
        
class Linear2DInterpolator():

    def __init__(self, xgrid, ygrid, x, y):
        self.x = x
        self.y = y
        self.xi = (x - xgrid.min())/(xgrid.max()-xgrid.min())*len(xgrid)
        self.yi = (y - ygrid.min())/(ygrid.max()-ygrid.min())*len(ygrid)
        self.valid_arrays = []
        
    # "Grows" the unmasked areas by one pixel
    @staticmethod
    def expandData(in_data):
        if not isinstance(in_data.mask, np.ndarray):
            in_data.mask = np.ones(in_data.shape, dtype=bool)*in_data.mask
        out_mask = ~ndimage.morphology.binary_dilation(~in_data.mask.copy())
        out_data = in_data.filled(np.finfo(np.float64).min)
        out_data = ndimage.morphology.grey_dilation(out_data, size=3)
        out_data[~in_data.mask] = in_data.data[~in_data.mask]
        out = np.ma.masked_array(out_data, mask=out_mask)
        return out

    def __call__(self, array2d):
        if not isinstance(array2d,np.ma.MaskedArray):
            logging.debug('Array used for interpolation is not a masked array.')
    
        array2d_ptr, read_only = array2d.data.__array_interface__['data'] 
    
        # change the array2d to fill masked parts 
        # i.e., land with no flow
        if array2d_ptr not in self.valid_arrays:
            # Fill eight cells of masked data
            for i in range(8):
                array2d = self.expandData(array2d)
            self.valid_arrays.append(array2d_ptr)
        
        interp = map_coordinates(array2d, [self.yi, self.xi],
                                 cval=np.nan, order=1)
        invalid = np.where(interp.min() < -1e+10 or
                           interp.max() > 1e+10)[0]
        if len(invalid) > 0:
            logging.warning('Invalid values returned by LinearNDFast!')
            # Uncomment to visualise data extrapolation
            #import matplotlib.pyplot as plt
            #old = array2d.copy()
            #plt.subplot(2,1,1)
            #plt.imshow(old)
            #plt.subplot(2,1,2)
            #plt.plot(self.yi, self.xi, '*w')
            #plt.imshow(array2d)
            #plt.show()

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
        self.zi = np.round(z_interpolator(z)).astype(np.int)
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
        # Indices corresponding to layers above and below
        interp_zi = z_interpolator(z)
        self.index_above = np.floor(interp_zi).astype(np.int)
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


###########################
# ReaderBlock
###########################

class ReaderBlock():
    """Class to store and interpolate the output from a reader."""

    def __init__(self, data_dict,
                 interpolation_horizontal='ndimage',
                 interpolation_vertical='linear'):

        # Make pointers to data values, for convenience
        self.x = data_dict['x']
        self.y = data_dict['y']
        self.time = data_dict['time']
        self.data_dict = data_dict
        del self.data_dict['x']
        del self.data_dict['y']
        del self.data_dict['time']
        try:
            self.z = data_dict['z']
            del self.data_dict['z']
        except:
            self.z = None

        # Mask any extremely large values, e.g. if missing netCDF _Fill_value
        for var in self.data_dict:
            if isinstance(self.data_dict[var], np.ma.core.MaskedArray):
                self.data_dict[var] = np.ma.masked_outside(
                    self.data_dict[var], -1E+9, 1E+9)

        # Set 1D (vertical) and 2D (horizontal) interpolators
        try:
            self.Interpolator2DClass = \
                horizontal_interpolation_methods[interpolation_horizontal]
        except Exception:
            raise NotImplementedError(
                'Valid interpolation methods are: ' +
                str(horizontal_interpolation_methods.keys()))

        try:
            self.Interpolator1DClass = \
                vertical_interpolation_methods[interpolation_vertical]
        except Exception:
            raise NotImplementedError(
                'Valid interpolation methods are: ' +
                str(vertical_interpolation_methods.keys()))

    def _initialize_interpolator(self, x, y, z=None):
        logging.debug('Initialising interpolator.')
        self.interpolator2d = self.Interpolator2DClass(self.x, self.y, x, y)
        if self.z is not None and len(np.atleast_1d(self.z)) > 1:
            self.interpolator1d = self.Interpolator1DClass(self.z, z)

    def interpolate(self, x, y, z=None, variables=None,
                    profiles=[], profiles_depth=None):

        self._initialize_interpolator(x, y, z)

        env_dict = {}
        if profiles is not []:
            profiles_dict = {'z': self.z}
        for varname, data in self.data_dict.iteritems():
            horizontal = self._interpolate_horizontal_layers(data)
            if profiles is not None and varname in profiles:
                profiles_dict[varname] = horizontal
            if horizontal.ndim > 1:
                env_dict[varname] = self.interpolator1d(horizontal)
            else:
                env_dict[varname] = horizontal
        return env_dict, profiles_dict

    def _interpolate_horizontal_layers(self, data):
        '''Interpolate all layers of 3d (or 2d) array.'''

        if data.ndim == 2:
            return self.interpolator2d(data)
        if data.ndim == 3:
            num_layers = data.shape[0]
            # Allocate output array
            result = np.ma.empty((num_layers, len(self.interpolator2d.x)))
            for layer in range(num_layers):
                result[layer, :] = self.interpolator2d(data[layer, :, :])
            return result

    def covers_positions(self, x, y, z=None):
        '''Check if given positions are covered by this reader block.'''

        indices = np.where((x >= self.x.min()) & (x <= self.x.max()) &
                           (y >= self.y.min()) & (y <= self.y.max()))[0]
        if len(indices) == len(x):
            return True
        else:
            return False
