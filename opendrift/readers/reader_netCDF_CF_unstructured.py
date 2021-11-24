# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2020, Gaute Hope, MET Norway
# Copyright 2020, Simon Weppe, MetOcean Solution, MetService New Zealand
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import numpy as np
from datetime import datetime, timedelta, timezone
from netCDF4 import Dataset, MFDataset, num2date
from scipy.interpolate import LinearNDInterpolator
import scipy, scipy.spatial
import pyproj
import logging
logger = logging.getLogger(__name__)

from opendrift.readers.basereader import BaseReader, UnstructuredReader


class Reader(BaseReader, UnstructuredReader):
    """
    A reader for unstructured (irregularily gridded) `CF compliant
    <https://cfconventions.org/>`_ netCDF files.

    Args:
        :param filename: A single netCDF file, or a pattern of files. The
                         netCDF file can also be an URL to an OPeNDAP server.
        :type filename: string, requiered.

        :param name: Name of reader
        :type name: string, optional

        :param proj4: PROJ.4 string describing projection of data.
        :type proj4: string, optional

    .. seealso::

        py:mod:`opendrift.readers.basereader.unstructured`.
    """

    # Misspelled standard names in some (Akvaplan-NIVA) FVCOM files
    variable_aliases = {
        'eastward_sea_water_velocity': 'x_sea_water_velocity',
        'Northward_sea_water_velocity': 'y_sea_water_velocity',
        'eastward wind': 'x_wind',
        'northward wind': 'y_wind',
        'ww': 'upward_sea_water_velocity',
    }

    node_variables = [
        'salinity',
        'temperature',
        'sea_floor_depth_below_sea_level',
        'sea_floor_depth_below_geoid',
        'sea_surface_height_above_geoid',
    ]

    face_variables = [
        'x_sea_water_velocity',
        'y_sea_water_velocity',
        'upward_sea_water_velocity',
    ]

    dataset = None

    # For in-memory caching of Sigma-coordinates and ocean depth
    siglay = None
    siglev = None
    siglay_center = None
    siglev_center = None
    ocean_depth_nele = None
    ocean_depth_node = None

    def __init__(self, filename=None, name=None, proj4=None):
        if filename is None:
            raise ValueError('Filename is missing')
        filestr = str(filename)
        if name is None:
            self.name = filestr
        else:
            self.name = name

        # xarray currently does not handle this type of grid:
        # https://github.com/pydata/xarray/issues/2233

        self.timer_start("open dataset")
        logger.info('Opening dataset: ' + filestr)
        if ('*' in filestr) or ('?' in filestr) or ('[' in filestr):
            logger.info('Opening files with MFDataset')
            self.dataset = MFDataset(filename)
        else:
            logger.info('Opening file with Dataset')
            self.dataset = Dataset(filename, 'r')

        if proj4 is not None:
            logger.info('Using custom projection: %s..' % proj4)
            self.proj4 = proj4
        else:
            self.proj4 = self.dataset.CoordinateProjection.strip()
            if self.proj4 == 'none': self.proj4 = None
            logger.info('Reading projection..: %s', self.proj4)
            assert self.proj4 is not None and len(
                self.proj4
            ) > 0, "No projection in data-file, please supply to reader"

        logger.info('Reading grid and coordinate variables..')
        assert self.dataset.CoordinateSystem == "Cartesian", "Only cartesian coordinate systems supported"

        self.x = self.dataset['x'][:]
        self.y = self.dataset['y'][:]
        self.xc = self.dataset['xc'][:]
        self.yc = self.dataset['yc'][:]

        # Times in FVCOM files are 'Modified Julian Day (MJD)', or fractional days since
        # 1858-11-17 00:00:00 UTC
        assert self.dataset['time'].time_zone == 'UTC'
        assert self.dataset['time'].units == 'days since 1858-11-17 00:00:00'
        assert self.dataset['time'].format == 'modified julian day (MJD)'
        ref_time = datetime(1858, 11, 17, 00, 00,
                            00)  # TODO: include , tzinfo=timezone.utc)
        self.times = np.array([
            ref_time + timedelta(days=d.item())
            for d in self.dataset['time'][:]
        ])
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        # time steps are not constant

        self.xmin = np.min(self.x)
        self.xmax = np.max(self.x)
        self.ymin = np.min(self.y)
        self.ymax = np.max(self.y)

        self.variable_mapping = {}
        for var_name in self.dataset.variables:
            # skipping coordinate variables
            if var_name in [
                    'x', 'y', 'time', 'lon', 'lat', 'lonc', 'latc', 'siglay',
                    'siglev', 'siglay_center', 'siglev_center'
            ]:
                continue

            # We only use sea_floor at nodes
            if var_name in ['h_center']:
                continue

            var = self.dataset[var_name]
            if 'standard_name' in var.ncattrs():
                std_name = var.getncattr('standard_name')
                std_name = self.variable_aliases.get(std_name, std_name)
                self.variable_mapping[std_name] = str(var_name)

        self.variables = list(self.variable_mapping.keys())

        # Run constructor of parent Reader class
        super().__init__()

        self.boundary = self._build_boundary_polygon_(self.x.compressed(),
                                                      self.y.compressed())

        self.timer_start("build index")
        logger.debug("building index of nodes..")
        self.nodes_idx = self._build_ckdtree_(self.x, self.y)

        logger.debug("building index of faces..")
        self.faces_idx = self._build_ckdtree_(self.xc, self.yc)

        self.timer_end("build index")

        self.timer_end("open dataset")

    def plot_mesh(self):
        """
        Plot the grid mesh. Does not automatically show the figure.
        """
        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter(self.x, self.y, marker='x', color='blue', label='nodes')
        plt.scatter(self.xc,
                    self.yc,
                    marker='o',
                    color='red',
                    label='centroids')

        x, y = getattr(self.boundary, 'context').exterior.xy
        plt.plot(x, y, color='green', label='boundary')

        plt.legend()
        plt.title('Unstructured grid: %s\n%s' % (self.name, self.proj))
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')

    def get_variables(self,
                      requested_variables,
                      time=None,
                      x=None,
                      y=None,
                      z=None):
        """
        FVCOM Grid:

        FVCOM uses 'triangular prisms' for gridding. Some variables are defined
        on the faces of the triangles, while others at the node.

        `x` and `y` holds the positions of the node, while `xc` and `yc` holds
        the positions on the centroids/faces. The centroids/faces are also
        known as 'zonal', or elements (presumably as in finite element).

        .. note::

            Currently this reader does not really interpolate. It looks up the
            closest point in time and space.

        Each element has a lookup-table of its surrounding elements, this list can be
        used when looking up elements for the interpolator of an arbitrary
        point on the grid. The same goes for the nodes.

        Let E be number of elements and N be number of nodes.

        Relevant lookup-tables:

        nbe:        [3 x E]  elements surround each element
        nbve:       [9 x N]  elements surrounding each node, minimum 3
        nbsn:       [11 x N] nodes surrounding each node

        Variables:

        Face:
        * u
        * v

        Node:
        * temperature
        * salinity

        """

        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        z = np.atleast_1d(z)
        if len(z) == 1:
            z = z[0] * np.ones(x.shape)

        logger.debug("Requested variabels: %s, lengths: %d, %d, %d" %
                     (requested_variables, len(x), len(y), len(z)))

        requested_variables, time, x, y, z, _outside = \
            self.check_arguments(requested_variables, time, x, y, z)

        nearest_time, _time_before, _time_after, indx_nearest, _indx_before, _indx_after = self.nearest_time(
            time)

        logger.debug("Nearest time: %s" % nearest_time)

        node_variables = [
            var for var in requested_variables if var in self.node_variables
        ]
        face_variables = [
            var for var in requested_variables if var in self.face_variables
        ]

        assert (len(node_variables) + len(face_variables)
                ) == len(requested_variables), "missing variables requested"

        variables = {}

        if node_variables:
            logger.debug("Interpolating node-variables..")

            nodes = self._nearest_node_(x, y)
            assert len(nodes) == len(x)

            for var in node_variables:
                dvar = self.variable_mapping.get(var)
                logger.debug("Interpolating: %s (%s)" % (var, dvar))
                dvar = self.dataset[dvar]

                # sigma ind depends on whether variable is defined on sigma layer og sigma level
                if 'siglev' in dvar.dimensions:
                    sigma_ind = self.__nearest_node_sigma__(dvar, nodes, z)
                    assert len(sigma_ind) == len(z)

                    # Reading the smallest block covering the actual data
                    block = dvar[indx_nearest,
                                 slice(sigma_ind.min(),
                                       sigma_ind.max() + 1),
                                 slice(nodes.min(),
                                       nodes.max() + 1)]
                    # Picking the nearest value
                    variables[var] = block[sigma_ind - sigma_ind.min(),
                                           nodes - nodes.min()]
                else:  # no depth dimension
                    variables[var] = dvar[slice(nodes.min(), nodes.max() + 1)]

        if face_variables:
            logger.debug("Interpolating face-variables..")

            fcs = self._nearest_face_(x, y)
            assert len(fcs) == len(x)

            for var in face_variables:
                dvar = self.variable_mapping.get(var)
                logger.debug("Interpolating: %s (%s)" % (var, dvar))
                dvar = self.dataset[dvar]

                # sigma ind depends on whether variable is defined on sigma layer og sigma level
                if 'siglay' in dvar.dimensions:
                    sigma_ind = self.__nearest_face_sigma__(dvar, fcs, z)
                    assert len(sigma_ind) == len(z)

                    # Reading the smallest profile slice covering the actual data
                    block = dvar[indx_nearest,
                                 slice(sigma_ind.min(),
                                       sigma_ind.max() + 1),
                                 slice(fcs.min(),
                                       fcs.max() + 1)]
                    # Picking the nearest value
                    variables[var] = block[sigma_ind - sigma_ind.min(),
                                           fcs - fcs.min()]
                else:  # no depth dimension
                    variables[var] = dvar[slice(fcs.min(), fcs.max() + 1)]

        return variables

    @staticmethod
    def _vector_nearest_(X, xp):
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
        return np.argmin(np.abs(X - xp), axis=0)

    def __nearest_node_sigma__(self, var, nodes, z):
        """
        Find nearest depth at node (sigma layer or sigma level depending on where the variable is defined).
        """
        shp = var.shape

        if self.siglay is None:
            logger.debug('Reading siglays into memory...')
            self.siglay = self.dataset['siglay'][:]

        if self.siglev is None:
            logger.debug('Reading siglevs into memory...')
            self.siglev = self.dataset['siglev'][:]

        if shp[1] == self.siglay.shape[0]:
            sigmas = self.siglay[:, nodes]
        else:
            sigmas = self.siglev[:, nodes]

        if self.ocean_depth_node is None:
            logger.debug('Reading ocean depth into memory...')
            self.ocean_depth_node = self.dataset['h'][:]

        # Calculating depths from sigmas
        depths = self.z_from_sigma(sigmas, self.ocean_depth_node[nodes])

        return self._vector_nearest_(depths, z)

    def __nearest_face_sigma__(self, var, el, z):
        """
        Find nearest depth at element (sigma layer or sigma level depending on where the variable is defined).
        """
        shp = var.shape

        if self.siglay_center is None:
            logger.info('Reading siglay_centers into memory...')
            self.siglay_center = self.dataset['siglay_center'][:]

        if self.siglev_center is None:
            logger.info('Reading siglev_centers into memory...')
            self.siglev_center = self.dataset['siglev_center'][:]

        if shp[1] == self.siglay_center.shape[0]:
            sigmas = self.siglay_center[:, el]
        else:
            sigmas = self.siglev_center[:, el]

        if self.ocean_depth_nele is None:
            logger.info('Reading ocean depth center into memory...')
            self.ocean_depth_nele = self.dataset['h_center'][:]

        # Calculating depths from sigmas
        depths = self.z_from_sigma(sigmas, self.ocean_depth_nele[el])

        return self._vector_nearest_(depths, z)

    @staticmethod
    def z_from_sigma(sigma, depth, elevation=0):
        """
        Calculate z-depth from sigma constant.

        https://rdrr.io/github/edwardlavender/fvcom.tbx/src/R/depth_from_known.R

        Args:

            sigma       Sigma coefficient(s)
            depth       Depth below mean sea-surface
            elevation   Elevation of sea-surface (e.g. tidal)

        Returns: z, depth below sea-surface in meters.
        """
        return sigma * (depth + elevation)

