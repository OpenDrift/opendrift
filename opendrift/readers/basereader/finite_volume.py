from abc import abstractmethod
import numpy as np
import logging

logger = logging.getLogger(__name__)

from .variables import Variables
from opendrift.readers.interpolation.finite_volume import FiniteVolumeReaderBlock
from numba import njit, prange

class FiniteVolumeReader(Variables):
    """
    An unstructured ReaderBlock for Finite Volume ocean models with irregularly gridded data on an Arakawa B-grid
    
    The initial type of grid that this class supports are `triangular prisms
    <https://en.wikipedia.org/wiki/Types_of_mesh#Triangular_prism>`_.
    Unstructured in xy-coordinates, x and y is constant in z. z might be
    non-cartesian (e.g. sigma-levels).

    .. seealso::

        :py:mod:`opendrift.readers`

        :class:`.structured.StructuredReader`
    """

    # Number of parallel threads to use when possible
    boundary = None

    # nodes
    x, y = None, None
    node_variables = None  # list of std-name variables defined at nodes

    # faces
    xc, yc = None, None
    face_variables = None  # list of std-name variables defined at center of faces

    @abstractmethod
    def get_variables(self, variables, time=None, x=None, y=None, z=None):
        """
        Obtain and return values of the requested variables at all positions
        (x, y, z) closest to given time.
        """

    def _get_variables_interpolated_(self, variables, profiles, profiles_depth, time, x, y, z):
        """
        The function called by OpenDrift to get an estimate of velocities at the particle at timestep n.
        - makes use of thredding to access the next timestep
        - treats static and dynamic variables differently
        """
        logger.debug(f"Requested variabels: {variables}")

        x, y, z = np.atleast_1d(x), np.atleast_1d(y), np.atleast_1d(z)

        # Isn't this a bit assumptuous?
        if len(z) == 1:
            z = z[0] * np.ones(x.shape)

        requested_variables, time, x, y, z, _outside = self.check_arguments(variables, time, x, y, z)

        # Figure out where we are in this netCDF files timeline
        # ----
        (nearest_time, time_before, time_after,
         indx_nearest, _indx_before, _indx_after,
        ) = self.nearest_time(time)
        logger.debug("Reader time:\n\t\t%s (before)\n\t\t%s (after)" % (time_before, time_after))
    

        # If none of the requested variables are dynamic, we won't need the blocks.
        # For time-dependent variables, however, we do.
        # ----
        time_dependent = [var for var in requested_variables if var not in self.static_variables]
        if not time_dependent and self.block_before is not None:
            env_before, env_after = self._interpolate_from_block_to_particles(
                x, y, z, time, time_before, time_after, variables
            )
            env = env_before
        else:
            # Get data at timesteps before- and after time (these are added to self)
            # ---
            logger.debug(f"interpolating {requested_variables} to the particles")
            self.timer_start("load_to_blocks")
            self._get_var_blocks(variables, time_before, time_after, x, y, z)
            self._check_var_blocks(variables, time_before, time_after, x, y, z)
            self.timer_end("load_to_blocks")
    
            # Interpolate data at those timesteps to the particles
            self.timer_start("interpolation_in_space")
            env_before, env_after = self._interpolate_from_block_to_particles(
                x, y, z, time, time_before, time_after, variables
            )
            self.timer_end("interpolation_in_space")
    
            # Interpolation in time
            self.timer_start("interpolation_in_time")
            if (
                (time_after is not None)
                and (time_before != time)
                and self.always_valid is False
            ):
                weight_after = self._find_weight_after(time)
                weight_before = 1-weight_after
                logger.debug(
                    ("Interpolating before (%s, weight %.2f) and\n\t\t      after (%s, weight %.2f) in time")
                    % (self.block_before.time, weight_before, self.block_after.time, weight_after)
                )
    
                # Interpolation procedure
                env = {}
                for var in variables:
                    env[var] = env_before[var]*weight_before + env_after[var]*weight_after
            else:
                logger.debug("No time interpolation needed - right on time.")
                env = env_before
            self.timer_end("interpolation_in_time")
        env_profiles = None
        return env, env_profiles

    def _find_weight_after(self, time):
        '''Calculates weight for interpolation in time'''
        weight_after = ((time - self.block_before.time).total_seconds()
                        / (self.block_after.time - self.block_before.time).total_seconds())

        assert weight_after <= 1, f'{weight_after=} is invalid'
        assert weight_after >= 0, f'{weight_after=} is invalid'
        return weight_after

    # blocks of data for timesteps
    # ----
    def _get_var_blocks(self, variables, time_before, time_after, x, y, z):
        """
        Make block variables for timestep before and after
            - variables    -- to be interpolated from the ocean model
            - time_before  -- model timestep before particle timestep
            - time_after   -- model timestep after particle timestep
            - x, y, z      -- Currently does not do anything in particular, since profiles are not supported
                              and interpolation is done later.

                              However get_variables (further down the pipeline) is expected to read x,y,z, so I just kept
                              it here to keep the call consistent with what OpenDrift expects upstream.
        """
        # Timestep (assuming times are regularly spaced)
        # ----
        dt = time_after-time_before

        # Swap block when moving over to the next model timestep, see if we can access the next block
        # ----
        if self.block_before is not None and self.block_after is not None:
            if self.block_before.time != time_before and self.block_after.time == time_before:
                self.block_before = self.block_after
                next_block = self.get_block_next(time_after)
                if next_block is not None:
                    self.block_after = next_block
                logger.debug(f"Fetched env-block for time before {self.block_before.time}")
                self.make_block_next(variables, time_after + dt, x, y, z)

            if self.block_after.time != time_after and self.block_before.time == time_after:
                self.block_after  = self.block_before
                next_block = self.get_block_next(time_before)
                if next_block is not None:
                    self.block_before = next_block
                logger.debug(f"Fetched env-block for time after {self.block_after.time}")
                self.make_block_next(variables, time_before - dt, x, y, z)

        # First access to blocks - when no buffer is available
        # ----
        if self.block_before is None or self.block_before.time != time_before:
            _ = self.get_block_next(time_before) # make sure to close the next_block thread if active
            self.block_before = self._fetch_data(variables, time_before, x, y, z)
            logger.debug(f"Fetched env-block for time before {self.block_before.time}")

        if self.block_after is None or self.block_after.time != time_after:
            _ = self.get_block_next(time_after) # make sure to close the next_block thread if active
            self.block_after = self._fetch_data(variables, time_after, x, y, z)
            logger.debug(f"Fetched env-block for time after {self.block_after.time}")

    def _check_var_blocks(self, variables, time_before, time_after, x, y, z):
        """
        check that both blocks have the data they need to update the environent
        """
        # Check that all requested variables are in the blocks, update both if something is missing
        if self.block_before is not None:
            missing = [var for var in variables if var not in self.block_before.data_dict.keys()]
            updated_variables = [var for var in self.block_before.data_dict.keys() if var != 'time'] + missing
            if missing:
                logger.info(f'Missing {missing}, reloading block_before with them')
                _ = self.get_block_next(time_before)
                self.block_before = self._fetch_data(updated_variables, time_before, x, y, z)
        
        if self.block_after is not None:
            missing = [var for var in variables if var not in self.block_after.data_dict.keys()]
            updated_variables = [var for var in self.block_before.data_dict.keys() if var != 'time'] + missing
            if missing:
                logger.info(f'Missing {missing}, reloading block_after with them')
                _ = self.get_block_next(time_after)
                self.block_after  = self._fetch_data(updated_variables, time_after,  x, y, z)
    
    def _fetch_data(self, variables, time, x, y, z):
        """
        Get another block of data if not already cached
        """
        var_block = FiniteVolumeReaderBlock(self.get_variables(variables, time, x, y, z))

        # Add some necessary variable info to the object
        var_block.node_variables   = self.node_variables
        var_block.face_variables   = self.face_variables
        var_block.variable_mapping = self.variable_mapping
        var_block.variable_aliases = self.variable_aliases
        return var_block

    # Methods needed to use an asynchronous strategy to speed up the reading part
    # ----
    def get_block_next(self, time):
        '''
        time is the timestep we want to read.
        Also used to ensure that the thread is reading data from a netCDF file before accessing more data in the netCDF
        '''
        block_next = None
        if self.block_next is not None:
            block_next = self.block_next.get()
            if block_next.time != time:
                block_next = None
            else:
                logger.debug(f"Fetched env-block for next time {block_next.time}")
        return block_next

    def make_block_next(self, variables, next_time, x, y, z):
        '''
        Start reading the (expected) next block of data
        '''
        if self.end_time > next_time:
            logger.debug(f"Fetcing env-block for next timestep {next_time}")
            _ = self.get_block_next(next_time) # just to make sure that all thredds are closed
            self.block_next = self.pool.apply_async(self._fetch_data, args = (variables, next_time, x, y, z))
        else:
            self.block_next = None
            
    # Interpolation
    # ----
    def _interpolate_from_block_to_particles(self, x, y, z, time, time_before, time_after, variables):
        """
        Interpolate data from block to positions of particles in space
        - Warning: Does not interpolate. The routine finds nearest point in space to (x, y, z).
        """
        node_variables = [var for var in variables if var in self.node_variables]
        cell_variables = [var for var in variables if var in self.face_variables]

        # Find nearest nodes, cells and corresponding sigma layers
        # ----
        cells = {}
        if cell_variables:
            cells["id"] = self._nearest_cell_(x, y) # 0.1 s for 1 millon particles, 2 million cells
            if self.use_3d:
                cells["sigma_ind"] = self.__nearest_cell_sigma__(cells["id"], z) # 0.4 s (same exp)
                (cells["sigma_next"], cells['weight_sigma_ind'], cells['weight_sigma_next'],
                    ) = self._sigma_for_interpolation(self.ocean_depth_cells, self.siglay_center, cells, z) # 0.2 s (same exp)

        nodes = {}
        if node_variables:
            nodes["id"] = self._nearest_node_(x, y) 
            if self.use_3d:
                nodes["sigma_ind"] = self.__nearest_node_sigma__(nodes["id"], z)
                (nodes["sigma_next"], nodes['weight_sigma_ind'], nodes['weight_sigma_next'],
                    ) = self._sigma_for_interpolation(self.ocean_depth_nodes, self.siglay, nodes, z)
        logger.debug(f"Interpolating before {self.block_before.time} in space {self.interpolation}")

        # Get environment at timestep before
        env_before = self.block_before.interpolate(nodes, cells, variables)

        # If relevant, get environment at timestep after
        if (time_after is not None) and (time_before != time):
            logger.debug(f"Interpolating after {self.block_after.time} in space {self.interpolation}")
            env_after = self.block_after.interpolate(nodes, cells, variables)
        else:
            env_after = env_before
        return env_before, env_after

    def __nearest_node_sigma__(self, nodes, z):
        """
        Find nearest depth at node.
        """
        # For now, we only need siglay
        sigmas = self.siglay[:, nodes]

        # Calculating depths from sigmas
        depths = self.z_from_sigma(sigmas, self.ocean_depth_nodes[nodes])

        return self._vector_nearest_(depths, z)

    def __nearest_cell_sigma__(self, cells, z):
        """
        Find nearest depth at element.
        """
        # add siglay if needed
        sigmas = self.siglay_center[:, cells]

        # Calculating depths from sigmas
        depths = self.z_from_sigma(sigmas, self.ocean_depth_cells[cells])

        return self._vector_nearest_(depths, z)

    def _sigma_for_interpolation(self, depth, sigma, grid, z):
        '''
        Find the sigma level which closes the "box" bounding nearest sigma and the particle
        - nearest --> sigma nearest
        - next    --> sigma next to nearest that creates a closed box covering the particle

        Using the vector for full depth is not very efficient, and can be rewritten
        - This routine has not been properly performance tested, I expect there to be room for rewrites...
        '''
        # Get grid depth and depth of the nearest sigma layer
        full_depths = self.z_from_sigma(sigma[:, grid['id']], depth[grid['id']])
        sigma_ind_depth = full_depths[grid['sigma_ind'], np.arange(grid['id'].shape[0])] # depth of nearest sigma level

        ## Determine if nearest sigma level is below or above the particle
        dz     = sigma_ind_depth-z

        # - Negative number if the particle is above nearest sigma, positive if particle is below.
        #   The other sigma level we wish to extract is therefore either above or below (+1 or -1)
        dsig   = self.sign(dz.data).astype(int)
        sigmas_next = grid['sigma_ind'] + dsig

        # Find levels where the particle is outside of the "sigma range"
        particles_over_shallowest_sigma = np.where(sigmas_next<0)[0]
        particles_below_deepest_sigma   = np.where(sigmas_next==full_depths.shape[0])[0]

        # We need these invalid sigma-layers to have valid indices when we broadcast, set to nearest sigma
        sigmas_next[particles_below_deepest_sigma] = grid['sigma_ind'][particles_below_deepest_sigma]

        # We use a linear interpolation scheme in the vertical and set next-sigma
        # weights to 0 where we just want to use the nearest sigma
        weight_sigma_next = np.abs(sigma_ind_depth - z)/np.abs(sigma_ind_depth - full_depths[sigmas_next, np.arange(grid['id'].shape[0])])
        weight_sigma_next[particles_below_deepest_sigma]   = 0
        weight_sigma_next[particles_over_shallowest_sigma] = 0

        assert weight_sigma_next.min() >= 0, f'invalid vertical interpolation weight: {weight_sigma_next.min()}'
        assert weight_sigma_next.max() <= 1, f'invalid vertical interpolation weight: {weight_sigma_next.max()}'

        weight_sigma_ind = 1-weight_sigma_next

        return sigmas_next, weight_sigma_ind, weight_sigma_next

    def _vector_nearest_(self, X, xp):
        """
        Find nearest element in vector of vectors `X` for each `xp`.

        Args:
            X   NxM matrix of levels
            xp  M   vector of positions

        Returns:
            i   M   vector of indices [0..N] of closest element in
                    X[0..N, i] to xp[i]
        """
        xp = np.atleast_2d(xp)
        diff = self._abs(X.data - xp)
        return self._find_sigma(diff)

    @staticmethod
    @njit
    def sign(field):
        return np.sign(field)

    @staticmethod
    @njit
    def _abs(difference):
        return np.abs(difference)

    @staticmethod
    @njit(parallel = True)
    def _find_sigma(depth_difference):
        """
        Essentially a sped-up version of np.argmin
        Benchmarked on 1 mill particles, 34 sigma layers:
            np.argmin(depth_difference, axis=0): 330 ms
            _find_sigma(depth_difference):        10 ms

        Developed and tested using numba-0.56.4
        """
        out = np.zeros((depth_difference.shape[1],), dtype=np.int32)
        for n in prange(depth_difference.shape[1]):
            out[n] = np.argmin(depth_difference[:, n])
        return out

    @staticmethod
    def z_from_sigma(sigma, depth, elevation=0):
        """
        Calculate z-depth from sigma constant.

        Args:
            sigma       Sigma coefficient(s)
            depth       Depth below mean sea-surface
            elevation   Elevation of sea-surface (e.g. tidal)

        Returns: z, depth below sea-surface in meters.
        """
        return sigma * (depth + elevation)

    def _nearest_node_(self, x, y):
        """
        Return nearest node (id) for x and y
        """
        query_points = np.array([x, y], dtype=np.float32).T
        _, inds = self.KDTree_node.query(query_points, k=1)
        return inds

    def _nearest_cell_(self, x, y):
        """
        Return nearest cell or face (id) for x and y
        """
        query_points = np.array([x, y], dtype=np.float32).T
        _, inds = self.KDTree_cell.query(query_points, k=1)
        return inds

    # Rotines expected by OpenDrift
    # ----
    def _build_boundary_polygon_(self, x_grid, y_grid):
        """
        Build a mesh boundary polygon
        """
        from shapely.geometry import Polygon
        from shapely.prepared import prep
        from scipy.spatial import ConvexHull

        P = np.vstack((x_grid, y_grid)).T
        hull = ConvexHull(P)

        boundary = P[hull.vertices, :]
        boundary = prep(Polygon(boundary))

        return boundary

    def covers_positions(self, x, y, z=0):
        """
        Check which points are within boundary of mesh.
        """
        assert (
            self.boundary is not None
        ), "Boundary of mesh has not been prepared by reader"

        # TODO: Check z coordinates
        logger.warning("z-coordinates are not bounds-checked")

        from shapely.vectorized import contains

        return contains(self.boundary, x, y)
