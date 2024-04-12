from abc import abstractmethod
import numpy as np
import scipy
import logging
logger = logging.getLogger(__name__)

from .variables import Variables


class UnstructuredReader(Variables):
    """
    An unstructured reader. Data is gridded irregularily.

    The initial type of grid that this class supports are `triangular prisms
    <https://en.wikipedia.org/wiki/Types_of_mesh#Triangular_prism>`_.
    Unstructured in xy-coordinates, x and y is constant in z. z might be
    non-cartesian (e.g. sigma-levels).

    .. seealso::

        :py:mod:`opendrift.readers`

        :class:`.structured.StructuredReader`

    """
    # Number of parallel threads to use when possible
    PARALLEL_WORKERS = -1

    boundary = None

    # nodes
    x = None
    y = None
    node_variables = None # list of std-name variables defined at nodes
    nodes_idx = None

    # faces
    xc = None
    yc = None
    face_variables = None # list of std-name variables defined at center of faces
    faces_idx = None

    def __init__(self):
        super().__init__()

    @abstractmethod
    def get_variables(self, variables, time=None, x=None, y=None, z=None):
        """
        Obtain and return values of the requested variables at all positions
        (x, y, z) closest to given time.

        Returns:

            Dictionary with arrays of length len(time) with values at exact positions x, y and z.
        """

    def _get_variables_interpolated_(self, variables, profiles,
                                   profiles_depth, time,
                                   reader_x, reader_y, z):

        env = self.get_variables(variables, time, reader_x, reader_y, z)

        # We probably have to use an UnstructuredBlock to store closest time-steps, and thus avoid fetching
        # more data on every call.

        logger.debug('Fetched env-before')
        env_profiles = None
        if profiles is not None:
            # Copying data from environment to vertical profiles
            env_profiles = {'z': [0, -profiles_depth]}
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

        .. note::

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

    def _build_rtree_(self, x, y):
        """
        Builds an R-tree of x, y
        """
        from rtree import index

        data = ((i, (xy[0], xy[1], xy[0], xy[1]), None) for i, xy in enumerate(zip(x, y)))
        idx = index.Index(data)

        return idx

    def _build_ckdtree_(self, x, y):
        from scipy.spatial import cKDTree
        P = np.vstack((x, y)).T
        return cKDTree(P)

    def __nearest_ckdtree__(self, idx, x, y):
        """
        Return index of nearest point in cKDTree
        """
        q = np.vstack((x, y)).T
        return idx.query(q, k = 1, workers = self.PARALLEL_WORKERS)[1]

    @staticmethod
    def __nearest_rtree__(idx, x, y):
        """
        Take array of points and get nearest point in rtree index.
        """
        return [idx.nearest((x, y, x, y), objects = False) for x, y in zip(x, y)]

    def _nearest_node_(self, x, y):
        """
        Return nearest node (id) for x and y
        """
        return self.__nearest_ckdtree__(self.nodes_idx, x, y)

    def _nearest_face_(self, xc, yc):
        """
        Return nearest element or face (id) for xc and yc
        """
        return self.__nearest_ckdtree__(self.faces_idx, xc, yc)
