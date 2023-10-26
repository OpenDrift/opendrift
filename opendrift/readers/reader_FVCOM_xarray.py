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

import numpy as np
from datetime import datetime, timedelta, timezone
from netCDF4 import Dataset
import xarray as xr
from scipy.interpolate import LinearNDInterpolator
import scipy, scipy.spatial
import pyproj
import logging
from opendrift.readers.basereader import BaseReader, UnstructuredReader

logger = logging.getLogger(__name__)

class Reader(BaseReader, UnstructuredReader):
    """
    A reader for unstructured (irregularily gridded) `CF compliant
    <https://cfconventions.org/>`_ netCDF files.

    Args:
        :param filename: A pattern of  netCDF files or a zarr archive. The
                         netCDF file can also be an URL to an OPeNDAP server.
        :type filename: string, required.

        :param filetype: 'netcdf' or 'zarr'.
        :type filetype: string, required.

        :param name: Name of reader
        :type name: string, optional

        :param proj4: PROJ.4 string describing projection of data.
        :type proj4: string, optional

    .. seealso::

        py:mod:`opendrift.readers.basereader.unstructured`.
    """

    # Misspelled standard names in some (Akvaplan-NIVA) FVCOM files
    variable_aliases = {
        # FVCOM name (standard_name, or long_name if necessary): OpenDrift standard
        'eastward_sea_water_velocity': 'x_sea_water_velocity',
        'Northward_sea_water_velocity': 'y_sea_water_velocity',
        'eastward wind': 'x_wind',
        'northward wind': 'y_wind',
        'ww': 'upward_sea_water_velocity',
        'sea_floor_depth_below_geoid': 'sea_floor_depth_below_sea_level',
    }

    node_variables = [
        'sea_floor_depth_below_sea_level',
        'sea_water_salinity',
        'sea_water_temperature',
        'sea_surface_height_above_geoid',
    ]

    face_variables = [
        'x_wind',
        'y_wind',
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

    def __init__(self, filename, proj4, filetype='zarr', cache=None, name=None):
        """
        filename: zarr store or netcdf file
        proj4: PROJ4 string
        filetype: "zarr" or "netcdf"
        cache: either None (no caching) or int (size of cache in bytes). 
            Only respected for zarr stores.
        name: name of Reader.
        """
        filestr = str(filename) # if Path was supplied
        if name is None:
            self.name = filestr
        else:
            self.name = name

        self.timer_start("open dataset")
        logger.info('Opening dataset: ' + filestr)
        if filetype=='zarr':
            from zarr.storage import FSStore, LRUStoreCache
            store_nocache = FSStore(filestr)
            if cache is not None:
                store = LRUStoreCache(store_nocache, max_size=cache) # e.g. 2**30 = 1 GiB cache
            else:
                store = store_nocache
            self.dataset = xr.open_dataset(store, engine='zarr', drop_variables=['time'])
        elif filetype=='netcdf':
            self.dataset = open_fvcom_files_as_xarray(filename)

        logger.info('Using custom projection: %s..' % proj4)
        self.proj4 = proj4

        logger.info('Reading grid and coordinate variables...')
        if self.dataset.CoordinateSystem == "Cartesian":
            self.x = self.dataset['x'][:]
            self.y = self.dataset['y'][:]
            self.xc = self.dataset['xc'][:]
            self.yc = self.dataset['yc'][:]
        elif self.dataset.CoordinateSystem == "GeoReferenced":
            self.x = self.dataset['lon'][:]
            self.y = self.dataset['lat'][:]
            self.xc = self.dataset['lonc'][:]
            self.yc = self.dataset['latc'][:]
        else:
            raise RuntimeError("This coordinate system is not supported")

        ref_time = datetime(1858, 11, 17, 00, 00,
                            00)  # TODO: include , tzinfo=timezone.utc)
        self.times = np.array([
            ref_time + timedelta(days=Itime.item()+Itime2.item()/86400e3)
            for Itime, Itime2 in zip(
                self.dataset['Itime'][:],
                self.dataset['Itime2'][:],
            )
        ])
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        # time steps are not constant

        self.xmin = np.min(self.x.values)
        self.xmax = np.max(self.x.values)
        self.ymin = np.min(self.y.values)
        self.ymax = np.max(self.y.values)

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

            attrs = self.dataset[var_name].attrs
            if 'standard_name' in attrs:
                std_name = attrs['standard_name']
                # TBD: check against requested_variables before adding
                logger.debug(f'Adding FVCOM std_name {std_name}')
                std_name = self.variable_aliases.get(std_name, std_name)
                logger.debug(f'Using OpenDrift alias: {std_name}')
                self.variable_mapping[std_name] = str(var_name)
            elif 'long_name' in attrs:
                std_name = attrs['long_name']
                logger.debug(f'Adding FVCOM long_name {std_name}')
                std_name = self.variable_aliases.get(std_name, std_name)
                logger.debug(f'Using OpenDrift alias: {std_name}')
                self.variable_mapping[std_name] = str(var_name)

        self.variables = list(self.variable_mapping.keys())

        # Run constructor of parent Reader class
        super().__init__()

        self.boundary = self._build_boundary_polygon_(self.x, self.y)

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

        all_avail = (len(node_variables) + len(face_variables)
                ) == len(requested_variables)
        if not all_avail:
            logger.debug(f"Node variables {node_variables} and face variables {face_variables}")
            logger.debug(f"Requested variables {requested_variables}")
            logger.debug(f"Variable mapping {self.variable_mapping}")
        assert all_avail, "missing variables requested"

        variables = {}

        x = xr.DataArray(x, dims=['particles'])
        y = xr.DataArray(y, dims=['particles'])
        z = xr.DataArray(z, dims=['particles'])

        if node_variables:
            logger.debug("Interpolating node-variables..")

            nodes = xr.DataArray(self._nearest_node_(x, y), dims=['particles'])
            
            assert len(nodes) == len(x)

            for var in node_variables:
                dvar = self.variable_mapping.get(var)
                logger.debug("Interpolating: %s (%s)" % (var, dvar))
                dvar = self.dataset[dvar]

                if any([d.startswith('sig') for d in dvar.dims]):
                    # super hacky!
                    sigma_ind = xr.DataArray(self.__nearest_node_sigma__(dvar, nodes.values, z.values), dims=['particles'])
                    assert len(sigma_ind) == len(z)
                else:
                    sigma_ind = slice(None)

                variables[var] = dvar.isel(
                    time=indx_nearest,
                    siglay_dim=sigma_ind,
                    siglev_dim=sigma_ind,
                    # hacky - relies on the fact that at most one of siglay/siglev
                    # is defined on any given array
                    node=nodes,
                    missing_dims='ignore',
                ).values

        if face_variables:
            logger.debug("Interpolating face-variables..")

            neles = xr.DataArray(self._nearest_face_(x, y), dims=['particles'])
            
            assert len(neles) == len(x)

            for var in face_variables:
                dvar = self.variable_mapping.get(var)
                logger.debug("Interpolating: %s (%s)" % (var, dvar))
                dvar = self.dataset[dvar]

                if any([d.startswith('sig') for d in dvar.dims]):
                    # super hacky!
                    sigma_ind = xr.DataArray(self.__nearest_face_sigma__(dvar, neles.values, z.values), dims=['particles'])
                    assert len(sigma_ind) == len(z)
                else:
                    sigma_ind = slice(None)

                variables[var] = dvar.isel(
                    time=indx_nearest,
                    siglay_dim=sigma_ind,
                    siglev_dim=sigma_ind,
                    # hacky - relies on the fact that at most one of siglay/siglev
                    # is defined on any given array
                    nele=neles,
                    missing_dims='ignore',
                ).values

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
        if 'siglay' in X.dims:
            xp = xr.DataArray(xp, dims=['node']) # add dummy dimension for easier xr magic
            inds = (X-xp).argmin(dim='siglay').values
        elif 'siglay_center' in X.dims:
            xp = xr.DataArray(xp, dims=['nele']) # add dummy dimension for easier xr magic
            inds = (X-xp).argmin(dim='siglay_center').values

        return inds

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

        if any([d.startswith('siglay') for d in var.dims]):
            sigmas = self.siglay[:]
        elif any([d.startswith('siglev') for d in var.dims]):
            sigmas = self.siglev[:]
        else:
            raise ValueError

        if self.ocean_depth_node is None:
            logger.debug('Reading ocean depth into memory...')
            self.ocean_depth_node = self.dataset['h'][:]

        # Calculating depths from sigmas
        depths = self.z_from_sigma(sigmas, self.ocean_depth_node[nodes])

        return self._vector_nearest_(depths, z)

    def __nearest_face_sigma__(self, var, el, z):
        """
        TBD: This is extremely hacky...

        Find nearest depth at element (sigma layer or sigma level depending on where the variable is defined).
        """
        shp = var.shape

        if self.siglay_center is None:
            logger.info('Reading siglay_centers into memory...')
            self.siglay_center = self.dataset['siglay_center'][:]

        if self.siglev_center is None:
            logger.info('Reading siglev_centers into memory...')
            self.siglev_center = self.dataset['siglev_center'][:]

        if any([d.startswith('siglay') for d in var.dims]):
            sigmas = self.siglay_center[:]
        elif any([d.startswith('siglev') for d in var.dims]):
            sigmas = self.siglev_center[:]
        else:
            raise ValueError

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

def open_fvcom_files_as_xarray(fname, load=False, **kwarg):
    """
    Works around xarray limitation to load FVCOM grid data
    """
    from glob import glob
    with Dataset(sorted(glob(fname))[0]) as f:
        ds = xr.open_mfdataset(fname, decode_times=False, drop_variables=['siglay', 'siglev'], data_vars='minimal', **kwarg)
        if 'siglay' in ds.dims:
            ds = _rename_sigma(ds, f, 'siglay')
        if 'siglev' in ds.dims:
            ds = _rename_sigma(ds, f, 'siglev')
    times = ds.Itime + ds.Itime2/86400e3
    ds['time'] = _timestamps_from_days_since_epoch(times.values)
    for dim, v in [
        ('time', 'Itime'),
        ('time', 'Itime2'),
        ('node', 'lat'),
        ('node', 'lon'),
        ('nele', 'lonc'),
        ('nele', 'latc'),
        ('node', 'x'),
        ('node', 'y'),
        ('nele', 'xc'),
        ('nele', 'yc'),
    ]:
        if v in ds.data_vars:
            ds = ds.assign_coords({v: (dim, ds[v].values)})
    return ds

def _timestamps_from_days_since_epoch(days):
    """
    Convert from iterable of days since epoch to a list of pandas Timestamps
    """
    import pandas as pd
    return [pd.Timestamp('1858-11-17 00:00:00') + pd.Timedelta(days=d) for d in days]

def _rename_sigma(ds, f, name):
    """
    name: siglay or siglev
    """
    dim_name = name+'_dim'
    ds = ds.rename_dims({name: dim_name})
    ds = ds.assign(variables={
        name: ((dim_name, 'node'), np.copy(f[name])),
    })
    ds = ds.assign_coords({
        name: ds[name].isel(node=0).values,
        name+'_center': ds[name+'_center'].isel(nele=0).values,
    })
    return ds
