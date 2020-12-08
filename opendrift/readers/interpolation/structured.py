import numpy as np
from scipy.ndimage import map_coordinates
import scipy.ndimage as ndimage
from scipy.interpolate import interp1d, LinearNDInterpolator

from .interpolators import Nearest2DInterpolator, fill_NaN_towards_seafloor, horizontal_interpolation_methods, vertical_interpolation_methods
from opendrift.readers.basereader import variables

import logging
logger = logging.getLogger(__name__)

class ReaderBlock():
    """Class to store and interpolate the output from a reader with data on a regular (structured) grid."""

    def __init__(self, data_dict,
                 interpolation_horizontal='linearNDFast',
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
        filled_variables = set()
        for var in self.data_dict:
            self.data_dict[var] = variables.Variables.__check_variable_array__(var, self.data_dict[var])

            if isinstance(self.data_dict[var], np.ma.core.MaskedArray):
                self.data_dict[var] = np.ma.masked_outside(
                    np.ma.masked_invalid(self.data_dict[var]), -1E+9, 1E+9)
                # Convert masked arrays to numpy arrays
                self.data_dict[var] = np.ma.filled(self.data_dict[var],
                                                   fill_value=np.nan)
            # Fill missing data towards seafloor if 3D
            if isinstance(self.data_dict[var], (list,)):
                logger.warning('Ensemble data currently not extrapolated towards seafloor')
            elif self.data_dict[var].ndim == 3:
                filled = fill_NaN_towards_seafloor(self.data_dict[var])
                if filled is True:
                    filled_variables.add(var)

        if len(filled_variables) > 0:
            logger.debug('Filled NaN-values toward seafloor for :'
                          + str(list(filled_variables)))

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

        if 'land_binary_mask' in self.data_dict.keys() and \
                interpolation_horizontal != 'nearest':
            logger.debug('Nearest interpolation will be used '
                          'for landmask, and %s for other variables'
                          % interpolation_horizontal)

    def _initialize_interpolator(self, x, y, z=None):
        logger.debug('Initialising interpolator.')
        self.interpolator2d = self.Interpolator2DClass(self.x, self.y, x, y)
        if self.z is not None and len(np.atleast_1d(self.z)) > 1:
            self.interpolator1d = self.Interpolator1DClass(self.z, z)

    def interpolate(self, x, y, z=None, variables=None,
                    profiles=[], profiles_depth=None):

        self._initialize_interpolator(x, y, z)

        env_dict = {}
        if profiles is not []:
            profiles_dict = {'z': self.z}
        for varname, data in self.data_dict.items():
            nearest = False
            if varname == 'land_binary_mask':
                nearest = True
                self.interpolator2d_nearest = Nearest2DInterpolator(self.x, self.y, x, y)
            if type(data) is list:
                num_ensembles = len(data)
                logger.debug('Interpolating %i ensembles for %s' % (num_ensembles, varname))
                if data[0].ndim == 2:
                    horizontal = np.zeros(x.shape)*np.nan
                else:
                    horizontal = np.zeros((len(self.z), len(x)))*np.nan
                ensemble_number = np.remainder(range(len(x)), num_ensembles)
                for en in range(num_ensembles):
                    elnum = ensemble_number == en
                    int_full = self._interpolate_horizontal_layers(data[en], nearest=nearest)
                    if int_full.ndim == 1:
                        horizontal[elnum] = int_full[elnum]
                    else:
                        horizontal[:, elnum] = int_full[:, elnum]
            else:
                horizontal = self._interpolate_horizontal_layers(data, nearest=nearest)
            if profiles is not None and varname in profiles:
                profiles_dict[varname] = horizontal
            if horizontal.ndim > 1:
                env_dict[varname] = self.interpolator1d(horizontal)
            else:
                env_dict[varname] = horizontal
        if 'z' in profiles_dict:
            profiles_dict['z'] = np.atleast_1d(profiles_dict['z'])

        return env_dict, profiles_dict

    def _interpolate_horizontal_layers(self, data, nearest=False):
        '''Interpolate all layers of 3d (or 2d) array.'''

        if nearest is True:
            interpolator2d = self.interpolator2d_nearest
        else:
            interpolator2d = self.interpolator2d
        if data.ndim == 2:
            return interpolator2d(data)
        if data.ndim == 3:
            num_layers = data.shape[0]
            # Allocate output array
            result = np.ma.empty((num_layers, len(interpolator2d.x)))
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

