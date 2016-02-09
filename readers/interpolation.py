import copy
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.interpolate import interp1d, LinearNDInterpolator


horizontal_interpolation_methods =  ['nearest', 'linearND', 'ndimage']
vertical_interpolation_methods =  ['nearest', 'linear']

class ReaderBlock():

    def __init__(self, data_dict, interpolation_method='ndimage'):

        # Make pointers to data values
        self.x = data_dict['x']
        self.y = data_dict['y']
        self.z = data_dict['z']
        self.time = data_dict['time']
        self.data_dict = copy.deepcopy(data_dict)
        del self.data_dict['x']
        del self.data_dict['y']
        del self.data_dict['z']
        del self.data_dict['time']

        if interpolation_method in horizontal_interpolation_methods:
            self.horizontal_interpolation_method = interpolation_method
        else:
            raise NotImplementedError('Valid interpolation methods '
                'are: ' + str(interpolation_methods))

    def interpolate(self, x, y, z=None, variables=None,
                    profiles=None, profiles_depth=None):

        #########################
        # Prepare interpolators
        #########################

        self.element_x = x
        self.element_y = y
        self.element_y = z

        # Nearest neighbour
        if self.horizontal_interpolation_method == 'nearest':
            xMin = self.x.min()
            xMax = self.x.max()
            yMin = self.y.min()
            yMax = self.y.max()
            self.xi = (x - xMin)/(xMax-xMin)*len(self.x)
            self.yi = (y - yMin)/(yMax-yMin)*len(self.y) 
            self.xi = np.round(self.xi).astype(np.int)
            self.yi = np.round(self.yi).astype(np.int)
            self.xi[self.xi>=len(self.x)] = len(self.x)-1
            self.yi[self.yi>=len(self.y)] = len(self.y)-1

            def interpolate_horizontal(slice2d):
                return slice2d[self.yi, self.xi]

        # scipy.ndimage
        if self.horizontal_interpolation_method == 'ndimage':
            xMin = self.x.min()
            xMax = self.x.max()
            yMin = self.y.min()
            yMax = self.y.max()
            self.xi = (x - xMin)/(xMax-xMin)*len(self.x)
            self.yi = (y - yMin)/(yMax-yMin)*len(self.y) 

            def interpolate_horizontal(slice2d):
                slice2d[slice2d.mask] = np.nan  # Gives holes
                return map_coordinates(slice2d, [self.yi, self.xi],
                                       cval=np.nan, order=0)

        # scipy.interpolate.LinearNDInterpolator
        elif self.horizontal_interpolation_method == 'linearND':
            block_x, block_y = np.meshgrid(self.x, self.y)
            self.block_x = block_x.ravel()
            self.block_y = block_y.ravel()
            self.interpolator = None

            def interpolate_horizontal(slice2d):
                valid = ~slice2d.ravel().mask

                if self.interpolator is None:
                    # Make new interpolator for given x,y
                    self.interpolator = LinearNDInterpolator(
                        (self.block_x[valid],
                         self.block_y[valid]),
                        slice2d.ravel()[valid])
                else:
                    # Reuse stored interpolator with new data
                    self.interpolator.values[:,0] = \
                        (slice2d.ravel()[valid])

                return self.interpolator(x, y)

        # Add the created interpolator function to self (Block)
        self.interpolate_horizontal = interpolate_horizontal
    
        if variables == 'prepare':
            return

        #############################
        # Interpolate all variables
        #############################
        env_dict = {}
        for varname, data in self.data_dict.iteritems():
            env_dict[varname] = \
                self.interpolate_horizontal_layers(data)
        return env_dict

    def interpolate_horizontal_layers(self, data):
        '''Interpolate all layers of 3d (or 2d) array.'''

        # Make a copy of the data array, to avoid surprices
        # due to multiple pointers
        data = copy.deepcopy(data)

        if data.ndim == 2:
            return np.ma.masked_invalid(self.interpolate_horizontal(data))
        if data.ndim == 3:
            num_layers = data.shape[0]
            # Allocate output array
            result = np.ma.empty((num_layers, len(self.element_x)))
            for layer in range(num_layers):
                result[layer,:] = \
                    self.interpolate_horizontal(data[layer,:,:])
            result = np.ma.masked_invalid(result)
        return result
