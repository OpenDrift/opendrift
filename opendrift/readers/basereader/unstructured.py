import numpy as np
import logging
logger = logging.getLogger(__name__)

from .variables import Variables

class UnstructuredReader(Variables):
    """
    An unstructured reader. Data is gridded irregularily.
    """

    var_block_before = None
    var_block_after  = None
    boundary = None

    def __init__(self):
        super().__init__()

        # Dictionaries to store blocks of data for reuse (buffering)
        self.var_block_before = {}  # Data for last timestep before present
        self.var_block_after = {}   # Data for first timestep after present

    def _get_variables_interpolated_(self, variables, profiles,
                                   profiles_depth, time,
                                   lon, lat, z,
                                   block, rotate_to_proj,
                                   ind_covered, reader_x, reader_y):

        env = self.get_variables_impl(variables, profiles,
                                            profiles_depth,
                                            time,
                                            reader_x, reader_y, z,
                                            block=block)

        logger.debug('Fetched env-before')
        env_profiles = None
        if profiles is not None:
            # Copying data from environment to vertical profiles
            env_profiles = {'z': profiles_depth}
            for var in profiles:
                env_profiles[var] = np.ma.array([env[var], env[var]])

        return env, env_profiles

    def _build_boundary_polygon_(self, x, y):
        """
        Build a polygon of the boundary of the mesh.

        Arguments:
            :param x: Array of node x position, lenght N
            :param y: Array of node y position, length N

        Returns:
            A `shapely.prepareped.prep` `shapely.Polygon`.

            The boundary of the mesh, ideally including holes in the mesh.

        Algorithms:
        -----------

        Try this alogrithm: https://stackoverflow.com/a/14109211/377927

        Boundary edges (line between two nodes) are only referenced by a single
        triangle.

        1. Find a starting edge segment: [v_start, v_next] (v is vertex or node)
        2. Find another _unvisited_ edge segment [v_i, v_j] that has
        either v_i = v_next or v_j = v_next and add the one not equal to v_next to the polygon.
        3. Reset v_next to the newly added point. Mark edge as visited.
        4. Continue untill we reach v_start.

        The polygon has a rotation, but this should not matter for our purpose
        of checking the bounds.

        Note: In order to find holes in the polygon all points must be scanned.

        Approximate using the convex hull:
        ----------------------
        An alternative simple approximation is to use the convex hull of the
        points, but this will miss points along the boundary which form a
        wedge in the boundary (as well as holes in the mesh).

        Holes in the mesh will often be covered by the landmask anyway, so they
        will usually not be a problem.

        """
        from shapely.geometry import Polygon
        from shapely.prepared import prep
        from scipy.spatial import ConvexHull

        P = np.vstack((x, y)).T
        hull = ConvexHull(P)

        boundary = P[hull.vertices, :]
        boundary = prep(Polygon(boundary))

        return boundary

    def covers_positions(self, x, y, z=0):
        """
        Check which points are within boundary of mesh.
        """
        assert self.boundary is not None, "Boundary of mesh has not been prepared by reader"

        # TODO: Check z coordinates
        logger.warning("z-coordinates are not bounds-checked")

        from shapely.vectorized import contains
        return contains(self.boundary, x, y)

