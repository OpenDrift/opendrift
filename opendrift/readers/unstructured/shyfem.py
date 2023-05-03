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
# Copyright 2021, Gaute Hope, MET Norway

import numpy as np
from datetime import datetime, timedelta
from netCDF4 import Dataset, MFDataset
import logging
logger = logging.getLogger(__name__)

from opendrift.readers.basereader import BaseReader, UnstructuredReader


class Reader(BaseReader, UnstructuredReader):
    """
    A reader for unstructured SHYFEM (irregularily gridded) `CF compliant
    <https://cfconventions.org/>`_ netCDF files.

    http://www.ismar.cnr.it/shyfem

    Args:
        :param filename: A single netCDF file, or a pattern of files. The
                         netCDF file can also be an URL to an OPeNDAP server.
        :type filename: string, requiered.

        :param name: Name of reader
        :type name: string, optional

    .. seealso::

        py:mod:`opendrift.readers.basereader.unstructured`.
    """

    variable_aliases = {
        'eastward_sea_water_velocity': 'x_sea_water_velocity',
        'northward_sea_water_velocity': 'y_sea_water_velocity',
        'sea_floor_depth_below_sea_surface': 'sea_floor_depth_below_sea_level'
    }

    dataset = None

    def __init__(self, filename=None, name=None):
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

        self.proj4 = '+proj=lonlat'

        logger.info('Reading grid and coordinate variables..')

        self.x, self.y = self.dataset['longitude'][:], self.dataset[
            'latitude'][:]

        ref_time = datetime.fromisoformat(self.dataset['time'].units[14:33])
        
        self.times = np.array([
            ref_time + timedelta(seconds=d.item())
            for d in self.dataset['time'][:]
        ])
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        # time steps are not constant

        self.xmin = np.min(self.x)
        self.xmax = np.max(self.x)
        self.ymin = np.min(self.y)
        self.ymax = np.max(self.y)

        # levels are the depth of the bottom of each layer. re-assign to middle of layer
        # for nearest interpolation.
        self.z = -self.dataset['level'][:]
        self.z = np.insert(self.z, 0, [0.])
        self.z = self.z[:-1] + (np.diff(self.z) / 2)
        assert len(self.z) == len(self.dataset['level'][:])
        self.zmin, self.zmax = np.min(self.z), 0.
        assert (self.z <= 0).all()

        self.variable_mapping = {}
        for var_name in self.dataset.variables:
            # skipping coordinate variables
            if var_name in ['time', 'longitude', 'latitude', 'levels']:
                continue

            var = self.dataset[var_name]
            if 'standard_name' in var.ncattrs():
                std_name = getattr(var, 'standard_name')
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
        self.timer_end("build index")

        self.timer_end("open dataset")

    def plot_mesh(self,corners=None):
        """
        Plot the grid mesh. Does not automatically show the figure.
        """
        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter(self.x, self.y, marker='x', color='blue', label='nodes')

        x, y = getattr(self.boundary, 'context').exterior.xy
        plt.plot(x, y, color='green', label='boundary')

        plt.legend()
        plt.title('Unstructured grid: %s\n%s' % (self.name, self.proj))
        plt.xlabel('lon [deg E]')
        plt.ylabel('lat [deg N]')

        if corners is not None:
            plt.xlim(corners[0],corners[1])
            plt.ylim(corners[2],corners[3])

    def get_variables(self,
                      requested_variables,
                      time=None,
                      x=None,
                      y=None,
                      z=None):
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        z = np.atleast_1d(z)
        #if len(z) == 1:
        #    z = z[0] * np.ones(x.shape)

        logger.debug("Requested variabels: %s, lengths: %d, %d, %d" %
                     (requested_variables, len(x), len(y), len(z)))

        requested_variables, time, x, y, z, _outside = \
            self.check_arguments(requested_variables, time, x, y, z)

        nearest_time, _time_before, _time_after, indx_nearest, _indx_before, _indx_after = self.nearest_time(
            time)

        logger.debug("Nearest time: %s" % nearest_time)

        variables = {}

        logger.debug("Interpolating node-variables..")

        nodes = self._nearest_node_(x, y)
        assert len(nodes) == len(x)

        for var in requested_variables:
            dvar = self.variable_mapping.get(var)
            logger.debug("Interpolating: %s (%s)" % (var, dvar))
            dvar = self.dataset[dvar]

            if len(dvar.shape) > 2:
                level_ind = self.__nearest_level__(z)

                # Reading the smallest block covering the actual data
                block = dvar[indx_nearest,
                             slice(nodes.min(),
                                   nodes.max() + 1),
                             slice(level_ind.min(),
                                   level_ind.max() + 1), ]

                # Picking the nearest value
                variables[var] = block[
                        nodes - nodes.min(),
                        level_ind - level_ind.min(),
                        ]
            elif len(dvar.shape) == 1:
                # Reading the smallest block covering the actual data
                block = dvar[slice(nodes.min(),
                                   nodes.max() + 1), ]

                # Picking the nearest value
                variables[var] = block[
                        nodes - nodes.min(),
                        ]
            else:
                logger.error('unknown dimensionality')

        return variables

    def __nearest_level__(self, z):
        """
        Find nearest index of z in levels.
        """
        return np.argmin(np.abs(self.z[:, None] - z), axis=0)
