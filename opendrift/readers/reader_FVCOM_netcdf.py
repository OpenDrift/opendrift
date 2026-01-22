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
# Copyright 2023, Håvard Espenes, University of Oslo / Akvaplan-niva
# Copyright 2020, Gaute Hope, MET Norway
# Copyright 2020, Simon Weppe, MetOcean Solution, MetService New Zealand
# Copyright 2015, Knut-Frode Dagestad, MET Norway

from datetime import datetime, timedelta
from netCDF4 import Dataset, MFDataset
from functools import cached_property
from pyproj import Proj
from pykdtree.kdtree import KDTree
from multiprocessing.pool import ThreadPool

import numpy as np
import logging

logger = logging.getLogger(__name__)

from opendrift.readers.basereader import BaseReader, FiniteVolumeReader

# This reader use threads to speed up reading from netCDF files. This, however, means that the
# reader can not access the dataset while OpenDrift is advecting particles (which would not be thread safe)
class Reader(BaseReader, FiniteVolumeReader):
    """
    A reader for unstructured (irregularily gridded) `CF compliant <https://cfconventions.org/>`_ netCDF files
    from models using the finite volume method. Note, we have just tested it with hydrodynamics from FVCOM 4.x.

    Args:
        :param filename: A single netCDF file, or a pattern of files. The netCDF file can also be an URL to an OPeNDAP server.
        :type filename: string, required.

        :param name: Name of reader
        :type name: string, optional

        :param proj4: PROJ.4 string describing projection of data.
        :type proj4: string, optional

    .. seealso::

        py:mod:`opendrift.readers.basereader.finite_volume`.
    """

    dataset = None
    block_before, block_after, block_next = None, None, None

    # We aspire to include more complicated schemes in the future.
    interpolation = "nearest"

    def __init__(
        self, filename=None, filelist=None, name=None, proj4=None, use_3d=True
    ):
        """
        The reader is used to read data, either from a single netCDF4 or from a number of files
        using the MFDataset class of netCDF4.

        - filename : use a specific file as input
        - filelist : use a list of files as input
        - name     : give a name to the reader
        - proj3    : projection of the readers grid coordinates
        - use_3d   : Switch
                     - if True (default): use (u,v,w) to advect particles
                     - if False:          use (ua, va) to advect particles
        """
        self.use_3d = use_3d
        assert (
            filename is not None or filelist is not None
        ), "You must either provide a filename or a filelist"
        filestr = str(filename)
        
        name
        if name is None:
            if filename is not None:
                self.name = filestr
            elif filelist is not None:
                self.name = filelist[0].split('/')[-1].split('_')[0]
        else:
            self.name = name
            
        self._open_fvcom_dataset(filestr, filename, filelist)

        self._load_fvcom_grid(proj4)

        self._check_and_prepare_time()

        self._configure_fields_to_read_from_fvcom()

        self._map_fvcom_variables_to_opendrift_internal_names()

        # Run constructor of parent Reader class
        super().__init__()

        self.boundary = self._build_boundary_polygon_(
            self.x.compressed(), self.y.compressed()
        )

        # Start a threadded pool for reading the next timestep
        self.pool = ThreadPool(processes = 1)

    @cached_property
    def KDTree_node(self):
        return KDTree(np.array([self.x, self.y], dtype=np.float32).T)

    @cached_property
    def KDTree_cell(self):
        return KDTree(np.array([self.xc, self.yc], dtype=np.float32).T)

    # Get bounds of the model domain
    # ----
    @property
    def xmin(self):
        return np.min(self.x)

    @property
    def xmax(self):
        return np.max(self.x)

    @property
    def ymin(self):
        return np.min(self.y)

    @property
    def ymax(self):
        return np.max(self.y)

    # Routines related to the initialization of the reader
    # ----
    def _open_fvcom_dataset(self, filestr, filename, filelist):
        """
        Opens fvcom file using netCDF dataset methods
        """
        self.timer_start("open dataset")

        if ("*" in filestr) or ("?" in filestr) or ("[" in filestr):
            logger.info(f"Opening multiple datasets")
            self.dataset = MFDataset(filename)

        elif filelist is not None:
            logger.info(f"Opening multiple datasets")
            self.dataset = MFDataset(filelist, "r")

        else:
            logger.info(f"Opening dataset: {filestr}")
            self.dataset = Dataset(filename, "r")
        self.timer_end("open dataset")

    def _load_fvcom_grid(self, proj4):
        """
        Makes the reader aware of projection and grid positions in this projection
        """
        logger.info("Reading grid and coordinate variables..")
        self.proj4 = proj4
        if self.dataset.CoordinateSystem == 'GeoReferenced':
            project = Proj(self.proj4)
            self.x, self.y   = project(self.dataset['lon'][:], self.dataset['lat'], inverse = True)
            self.xc, self.yc = project(self.dataset['lonc'][:], self.dataset['latc'], inverse = True)

        elif self.dataset.CoordinateSystem == "Cartesian":
            self.x, self.y   = self.dataset['x'][:], self.dataset['y'][:]
            self.xc, self.yc = self.dataset['xc'][:], self.dataset['yc'][:]

        # depth properties
        self.siglay, self.siglev = self.dataset['siglay'][:], self.dataset['siglev'][:]
        self.ocean_depth_nodes = self.dataset["h"][:]
        self.ocean_depth_cells = self.dataset["h_center"][:]
        self.siglay_center = self.dataset["siglay_center"][:]
        self.siglev_center = self.dataset["siglev_center"][:]

    def _check_and_prepare_time(self):
        """
        Checks that the time is in assumed format, defines the datum and updates the reader
        with min- and max time for the current file. This is essential information for the blocks later on.
        """
        assert self.dataset["time"].time_zone == "UTC"
        assert self.dataset["time"].units == "days since 1858-11-17 00:00:00"
        assert self.dataset["time"].format == "modified julian day (MJD)"
        ref_time = datetime(
            1858, 11, 17, 00, 00, 00
        )  # TODO: include , tzinfo=timezone.utc)
        self.times = np.array(
            [ref_time + timedelta(days=Itime.item()) + timedelta(milliseconds=Itime2.item()) for (Itime, Itime2) in zip(self.dataset["Itime"][:], self.dataset['Itime2'][:])]
        )
        self.start_time = self.times[0]
        self.end_time = self.times[-1]

    def _configure_fields_to_read_from_fvcom(self):
        """
        FVCOM stores fields with different long_names and variable names than done internally in
        opendrift. Furthermore, we want to be able to advect using either 2d or 3d fields.
        To do so, we need some dictionaries and lists that correctly references fields.

        Compatible with uk-fvcom branch (not sure how / if this differs from original FVCOM)
        """
        # A list of time-independent variables
        self.static_variables = [
            'sea_floor_depth_below_sea_level',
            'sea_surface_height_above_geoid',
        ]

        # And then make a list of all fields we might want to read from FVCOM
        if self.use_3d:
            logger.info("Load 3D (u, v, w) velocity data")
            # name_used_in_fvcom : equivalent_CF_name used in opendrift
            self.variable_aliases = {
                "eastward_sea_water_velocity": "x_sea_water_velocity",
                "Northward_sea_water_velocity": "y_sea_water_velocity",
                "eastward wind": "x_wind",
                "northward wind": "y_wind",
                "Upward Water Velocity": "upward_sea_water_velocity",
                "sea_floor_depth_below_geoid": "sea_floor_depth_below_sea_level",
                "bed stress magnitude from currents": "bed_stress",
            }

            self.node_variables = [
                "sea_water_salinity",
                "sea_water_temperature",
                "sea_floor_depth_below_sea_level",
                "sea_surface_height_above_geoid",
            ]

            self.face_variables = [
                "x_wind",
                "y_wind",
                "x_sea_water_velocity",
                "y_sea_water_velocity",
                "upward_sea_water_velocity",
            ]

            assert (
                "u" in self.dataset.variables.keys()
            ), f"u not found in {self.dataset.filepath()}"

            assert (
                "v" in self.dataset.variables.keys()
            ), f"v not found in {self.dataset.filepath()}"

        else:
            logger.info("Load 2D (ua, va) velocity data")
            # [name_used_in_fvcom : equivalent_CF_name_for_opendrift]
            self.variable_aliases = {
                "Vertically Averaged x-velocity": "x_sea_water_velocity",
                "Vertically Averaged y-velocity": "y_sea_water_velocity",
                "sea_floor_depth_below_geoid": "sea_floor_depth_below_sea_level",
                "bed stress magnitude from currents": "q",
            }

            self.node_variables = [
                "sea_surface_height_above_geoid",
                "sea_floor_depth_below_sea_level",
            ]

            self.face_variables = ["x_sea_water_velocity", "y_sea_water_velocity"]

            assert (
                "ua" in self.dataset.variables.keys()
            ), f"ua not found in {self.dataset.filepath()}"

            assert (
                "va" in self.dataset.variables.keys()
            ), f"va not found in {self.dataset.filepath()}"
    
    def _map_fvcom_variables_to_opendrift_internal_names(self):
        """
        since OpenDrift is a generic tool
        """
        self.variable_mapping = {}
        for var_name in self.dataset.variables:
            # skipping coordinate variables
            if var_name in [
                "x",
                "y",
                "time",
                "lon",
                "lat",
                "lonc",
                "latc",
                "siglay",
                "siglev",
                "siglay_center",
                "siglev_center",
            ]:
                continue

            var = self.dataset[var_name]
            if "standard_name" in var.ncattrs():
                std_name = var.standard_name
                std_name = self.variable_aliases.get(std_name, std_name)
                self.variable_mapping[std_name] = str(var_name)

            elif "long_name" in var.ncattrs():
                long_name = var.long_name
                long_name = self.variable_aliases.get(long_name, long_name)
                self.variable_mapping[long_name] = str(var_name)
        self.variables = list(self.variable_mapping.keys())

    def get_variables(self, requested_variables, time=None, x=None, y=None, z=None):
        """
        Routine used to get data from FVCOM results. FVCOM chunks by time, so it makes sense to extract all of it
        and dump to a data cacher (ReaderBlock)

        FVCOM Grid:
            - Uses triangular prisms

        Variables:
            - Face (vector variables):
            - Node (scalar variables)
        """
        logger.debug("Requested variabels: %s" % (requested_variables))
        requested_variables, time, x, y, z, _outside = self.check_arguments(
            requested_variables, time, x, y, z
        )
        (
            nearest_time,
            _time_before,
            _time_after,
            indx_nearest,
            _indx_before,
            _indx_after,
        ) = self.nearest_time(time)

        logger.debug("Nearest time: %s" % nearest_time)
        node_variables = [
            var for var in requested_variables if var in self.node_variables
        ]
        face_variables = [
            var for var in requested_variables if var in self.face_variables
        ]
        assert (len(node_variables) + len(face_variables)) == len(
            requested_variables
        ), "missing variables requested"

        # Let the routine using this dict know which time was read from FVCOM
        # ----
        variables = {"time": nearest_time}

        # Load data from FVCOM
        # ----
        if node_variables:
            logger.debug("Reading node-variables..")
            for var in node_variables:
                dvar = self.variable_mapping.get(var)
                dvar = self.dataset[dvar]
                logger.debug("Reading from netCDF: %s (%s)" % (var, dvar))
                if var == "nodes":
                    variables[var] = np.empty([])  # Hvorfor prøver den å lese nodes?
                elif dvar.ndim == 1:
                    variables[var] = dvar[:]
                else:
                    variables[var] = dvar[indx_nearest, :]

        if face_variables:
            logger.debug("Reading face-variables..")
            for var in face_variables:
                dvar = self.variable_mapping.get(var)
                logger.debug("Reading from netCDF: %s (%s)" % (var, dvar))
                dvar = self.dataset[dvar]
                if dvar.ndim == 1:
                    variables[var] = dvar[:]
                else:
                    variables[var] = dvar[indx_nearest, :]
        return variables
