import numpy as np
from abc import abstractmethod

from opendrift.timer import Timeable

from .consts import standard_names

class Variables(Timeable):
    derived_variables = None
    name = None

    def __init__(self):
        self.derived_variables = {}

    def __get_variables__(self, variables, profiles, profiles_depth,
                       time, x, y, z, block):
        """Wrapper around reader-specific function get_variables()

        Performs some common operations which should not be duplicated:
        - monitor time spent by this reader
        - convert any numpy arrays to masked arrays
        """
        self.logger.debug('Fetching variables from ' + self.name)
        self.timer_start('reading')

        if profiles is not None and block is True:
            # If profiles are requested for any parameters, we
            # add two fake points at the end of array to make sure that the
            # requested block has the depth range required for profiles
            x = np.append(x, [x[-1], x[-1]])
            y = np.append(y, [y[-1], y[-1]])
            z = np.append(z, [profiles_depth[0], profiles_depth[1]])
        env = self.get_variables_derived(variables, time, x, y, z, block)

        # Make sure x and y are floats (and not e.g. int64)
        if 'x' in env.keys():
            env['x'] = np.array(env['x'], dtype=np.float)
            env['y'] = np.array(env['y'], dtype=np.float)

        for variable in variables:
            # Convert any masked arrays to NumPy arrays
            if isinstance(env[variable], np.ma.MaskedArray):
                env[variable] = env[variable].filled(np.nan)
            # Mask values outside valid_min, valid_max (self.standard_names)
            if variable in standard_names.keys():
                if isinstance(env[variable], list):
                    self.logger.warning('Skipping min-max checking for ensemble data')
                    continue
                with np.errstate(invalid='ignore'):
                    invalid_indices = np.logical_and(
                        np.isfinite(env[variable]), np.logical_or(
                        env[variable]<standard_names[variable]['valid_min'],
                        env[variable]>standard_names[variable]['valid_max']))
                if np.sum(invalid_indices) > 0:
                    invalid_values = env[variable][invalid_indices]
                    self.logger.warning('Invalid values (%s to %s) found for %s, replacing with NaN' % (invalid_values.min(), invalid_values.max(), variable))
                    self.logger.warning('(allowed range: [%s, %s])' %
                                    (standard_names[variable]['valid_min'],
                                     standard_names[variable]['valid_max']))
                    env[variable][invalid_indices] = np.nan

        # Convolve arrays with a kernel, if reader.convolve is set
        if self.convolve is not None:
            from scipy import ndimage
            N = self.convolve
            if isinstance(N, (int, np.integer)):
                kernel = np.ones((N, N))
                kernel = kernel/kernel.sum()
            else:
                kernel = N
            self.logger.debug('Convolving variables with kernel: %s' % kernel)
            for variable in env.keys():
                if variable in ['x', 'y', 'z', 'time']:
                    pass
                else:
                    if env[variable].ndim == 2:
                        env[variable] = ndimage.convolve(
                            env[variable], kernel, mode='nearest')
                    elif env[variable].ndim == 3:
                        env[variable] = ndimage.convolve(
                            env[variable], kernel[:,:,None],
                            mode='nearest')

        self.timer_end('reading')

        return env

    def get_variables_derived(self, variables, *args, **kwargs):
        """Wrapper around get_variables, adding derived"""
        if isinstance(variables, str):
            variables = [variables]
        if not isinstance(variables, list):
            variables = list(variables)
        derive_variables = False
        for var in variables:
            if var in self.derived_variables:
                fromvars = self.derived_variables[var]
                for v in fromvars:
                    variables.append(v)
                # Removing the derived variable name
                variables = [v for v in variables if v != var]
                derive_variables = True

        env = self.get_variables(variables, *args, **kwargs)
        if derive_variables is True:
            self.calculate_derived_environment_variables(env)

        return env

    @abstractmethod
    def get_variables(self, variables, time=None,
                      x=None, y=None, z=None, block=False):
        """Method which must be implemented by all reader-subclasses.

        Obtain and return values of the requested variables at all positions
        (x, y, z) closest to given time.

        Arguments:
            variables: string, or list of strings (standard_name) of
                requested variables. These must be provided by reader.
            time: datetime or None, time at which data are requested.
                Can be None (default) if reader/variable has no time
                dimension (e.g. climatology or landmask).
            x, y: float or ndarrays; coordinates of requested points in the
                Spatial Reference System (SRS) _of the reader (NB!!)_.
            z: float or ndarray; vertical position (in meters, positive up)
                of requested points.
                default: 0 m (unless otherwise documented by reader)
            block: bool, see return below

          Returns:
            data: Dictionary
                keywords: variables (string)
                values:
                    - 1D ndarray of len(x) if block=False. Nearest values
                        (neighbour) of requested position are returned.
                    - 3D ndarray encompassing all requested points in
                        x,y,z domain if block=True. It is task of invoking
                        application (OpenDriftSimulation) to perform
                        interpolation in space and time.
        """

    @abstractmethod
    def get_variables_interpolated(self, variables, profiles=None,
                                   profiles_depth=None, time=None,
                                   lon=None, lat=None, z=None,
                                   block=False, rotate_to_proj=None):
        """
        `get_variables_interpolated` is the interface to the basemodel, and is
        responsible for returning variables at the correct positions. This is done by:

            1. Calling `__get_variables__` which,
            2. calls `get_variables_derived`, which, finally
            3. calls `get_variables`.

        This function is responsible for setting up an interpolator,
        `ReaderBlock` for regular gridded datasets.

        Arguments:
            variables: string, or list of strings (standard_name) of
                requested variables. These must be provided by reader.

            profiles: N/A

            profiles_depth: N/A

            time: datetime or None, time at which data are requested.
                Can be None (default) if reader/variable has no time
                dimension (e.g. climatology or landmask).

            lon: N/A

            lat: N/A

            z: float or ndarray; vertical position (in meters, positive up)
                of requested points.
                default: 0 m (unless otherwise documented by reader)

            block: bool, see return below

            rotate_to_proj: N/A

          Returns:
            data: Dictionary
                keywords: variables (string)
                values:
                    - 1D ndarray of len(x) if block=False. Nearest values
                        (neighbour) of requested position are returned.
                    - 3D ndarray encompassing all requested points in
                        x,y,z domain if block=True. It is task of invoking
                        application (OpenDriftSimulation) to perform
                        interpolation in space and time.
        """
        pass
