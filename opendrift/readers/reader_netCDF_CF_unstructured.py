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
# Copyright 2015, Knut-Frode Dagestad, MET Norway
# Copyright 2020, Simon Weppe, MetOcean Solution, MetService New Zealand
# Copyright 2020, Gaute Hope, MET Norway

import numpy as np
from datetime import datetime, timedelta, timezone
from netCDF4 import Dataset, MFDataset, num2date
from scipy.interpolate import LinearNDInterpolator
import pyproj
import logging
logger = logging.getLogger(__name__)

from opendrift.readers.basereader import BaseReader, UnstructuredReader

class Reader(BaseReader, UnstructuredReader):
    """
    A reader for unstructured (irregularily gridded) `CF-compliant
    <https://cfconventions.org/>`_ netCDF files.

    The reader interpolates on the projected grid (meters).

    Args:
        :param filename: A single netCDF file, or a pattern of files. The
                         netCDF file can also be an URL to an OPeNDAP server.
        :type filename: string, requiered.

        :param name: Name of reader
        :type name: string, optional

        :param proj4: PROJ.4 string describing projection of data.
        :type proj4: string, optional
    """

    # Misspelled standard names in some (Akvaplan-NIVA) FVCOM files
    variable_aliases = {
        'eastward_sea_water_velocity': 'x_sea_water_velocity',
        'Northward_sea_water_velocity': 'y_sea_water_velocity',
        'eastward wind': 'x_wind',
        'northward wind': 'y_wind'
        }

    # Mapping FVCOM variable names to CF standard_name
    fvcom_mapping = {
        'um': 'x_sea_water_velocity',
        'vm': 'y_sea_water_velocity'}

    dataset = None
    x = None
    y = None

    def __init__(self, filename = None, name = None, proj4 = None):
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
            assert self.proj4 is not None and len(self.proj4) > 0, "No projection in data-file, please supply to reader"

        logger.info('Reading grid and coordinate variables..')
        assert self.dataset.CoordinateSystem == "Cartesian", "Only cartesian coordinate systems supported"

        self.x = self.dataset['x'][:]
        self.y = self.dataset['y'][:]

        # Times in FVCOM files are 'Modified Julian Day (MJD)', or fractional days since
        # 1858-11-17 00:00:00 UTC
        assert self.dataset['time'].time_zone == 'UTC'
        assert self.dataset['time'].units == 'days since 1858-11-17 00:00:00'
        assert self.dataset['time'].format == 'modified julian day (MJD)'
        ref_time = datetime(1858, 11, 17, 00, 00, 00, tzinfo = timezone.utc)
        self.times = np.array([ ref_time + timedelta(days = d.item()) for d in self.dataset['time'][:] ])
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        # time steps are not constant

        self.xmin = np.min(self.x)
        self.xmax = np.max(self.x)
        self.ymin = np.min(self.y)
        self.ymax = np.max(self.y)

        self.variable_mapping = {}
        for var_name in self.dataset.variables:
            logger.debug("Parsing variable: %s" % var_name)

            # skipping coordinate variables
            if var_name in ['x', 'y', 'time', 'lon', 'lat', 'lonc', 'latc',
                    'siglay', 'siglev', 'siglay_center', 'siglev_center']:
                continue

            var = self.dataset[var_name]
            if 'standard_name' in var.ncattrs():
                std_name = var.getncattr('standard_name')
                std_name = self.variable_aliases.get(std_name, std_name)
                self.variable_mapping[std_name] = str(var_name)
            elif var_name in self.fvcom_mapping:
                std_name = self.fvcom_mapping[var_name]
                self.variable_mapping[std_name] = str(var_name)

        self.variables = list(self.variable_mapping.keys())

        # Run constructor of parent Reader class
        super().__init__()

        self.boundary = self._build_boundary_polygon_(self.x.compressed(), self.y.compressed())

        self.timer_end("open dataset")

    def plot_nodes(self):
        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter(self.x, self.y, marker = 'x', color = 'blue', label = 'nodes')
        plt.scatter(self.dataset['xc'][:], self.dataset['yc'][:], marker = 'o', color = 'red', label = 'centroids')

        x,y = getattr(self.boundary, 'context').exterior.xy
        plt.plot(x, y, color = 'green', label = 'boundary')

        plt.legend()
        plt.title('Unstructured grid: %s\n%s' % (self.name, self.proj))
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None, block=False):
        """
        FVCOM
        ---------------

        Grid
        ====

        FVCOM uses 'triangular prisms' for gridding. Some variables are defined
        on the faces of the triangles, while others at the node.

        `x` and `y` holds the positions of the node, while `xc` and `yc` holds
        the positions on the centroids/faces. The centroids/faces are also
        known as 'zonal', or elements (presumably as in finite element).


        Each element has a lookup-table of its surrounding elements, this list can be
        used when looking up elements for the interpolator of an arbitrary
        point on the grid. The same goes for the nodes.

        Let E be number of elements and N be number of nodes.

        Relevant lookup-tables
        =============

        cbe:        [3 x E]  elements surround each element
        nbve:       [9 x N]  elements surrounding each node, minimum 3
        nbsn:       [11 x N] nodes surrounding each node

        Variables
        =========

        Face:
        * u
        * v

        Node:
        * temperature
        * salinity

        """

        requested_variables, time, x, y, z, outside = \
            self.check_arguments(requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)


        x = np.atleast_1d(x)
        y = np.atleast_1d(y)


        pass
