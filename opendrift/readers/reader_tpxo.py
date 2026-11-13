"""
Reader for barotropic tidal currents from the TPXO9 global tide model, e.g.
to add tides to an ocean model without tides (such as TOPAZ5) with reader
arithmetic::

    tides = reader_tpxo.Reader('/path/to/tpxo9', corners=[-15, 40, 55, 86])
    o.add_reader([topaz + tides, topaz])

``topaz + tides`` provides only the variables both readers have (the
currents), so ``topaz`` is added again for its other variables. TPXO currents
are depth-averaged, so the same tidal current is returned at all depths.

The harmonic synthesis follows the OTIS/TMD convention TPXO is distributed for
(as pyTMD with corrections='OTIS'): constituent phase is
omega * (t - 1992-01-01) + phase0 + u, with nodal factors f and angles u from
the lunar node only, computed with pyTMD (which must be installed).

Reads the OTIS netCDF velocity file of TPXO9.v1 (``u_tpxo9.v1.nc``: ``ua``/``up``
eastward and ``va``/``vp`` northward velocity amplitude (cm/s) and Greenwich
phase lag (degrees) on a regular global lon/lat C-grid).
"""

import os
import logging
from datetime import datetime

import numpy as np

from opendrift.readers.basereader.continuous import ContinuousReader

logger = logging.getLogger(__name__)

TPXO_UV_FILE = 'u_tpxo9.v1.nc'

# The eight major semidiurnal and diurnal constituents
MAJOR_CONSTITUENTS = ('m2', 's2', 'n2', 'k2', 'k1', 'o1', 'p1', 'q1')

OTIS_EPOCH = datetime(1992, 1, 1)
MJD_OTIS_EPOCH = 48622.0


def tidal_phasors(time, constituents):
    """f * exp(i * (omega t + phase0 + u)) per constituent at ``time``
    (naive UTC datetime), as in pyTMD.predict.ocean_load.time_series with
    corrections='OTIS'. Multiply by a complex harmonic constant
    (amp * exp(-i * phase_lag)) and take the real part to get the tide."""
    import pyTMD
    days = (time - OTIS_EPOCH).total_seconds() / 86400.0
    u, f, _ = pyTMD.constituents.arguments(
        np.atleast_1d(days + MJD_OTIS_EPOCH), constituents, corrections='OTIS')
    _, phase0, omega, _, _ = pyTMD.constituents._constituent_parameters(constituents)
    return f[0] * np.exp(1j * (86400.0 * omega * days + phase0 + u[0]))


class HarmonicGrid:
    """Complex harmonic constants of one TPXO C-grid component on a regular
    lon/lat grid, bilinearly interpolated to arbitrary positions. Complex
    values (not amplitude/phase) are interpolated, so phase wrap-around
    doesn't matter; land nodes are 0 in TPXO and stay 0."""

    def __init__(self, ds, amp_var, phase_var, lon_var, lat_var, con_index,
                 corners=(-180, 180, -90, 90), amp_scale=0.01):
        # The file is (nc, nx, ny) with ny fastest: read whole longitude
        # blocks and cut the latitudes in memory, as reading a latitude band
        # directly means thousands of small strided reads (minutes on GPFS)
        lon_min, lon_max, lat_min, lat_max = corners
        lon = ds[lon_var][:][:, 0].astype(np.float64)
        lat = ds[lat_var][0, :].astype(np.float64)
        self.dlon = (lon[-1] - lon[0]) / (len(lon) - 1)
        self.dlat = (lat[-1] - lat[0]) / (len(lat) - 1)
        nlon = len(lon)
        j0 = max(int(np.searchsorted(lat, lat_min)) - 1, 0)
        j1 = min(int(np.searchsorted(lat, lat_max)) + 1, len(lat))
        nx = int(np.ceil((lon_max - lon_min) / self.dlon)) + 3
        self.periodic = nx >= nlon
        if self.periodic:
            i0, nx = 0, nlon
        else:
            i0 = (int(np.floor((lon_min - lon[0]) / self.dlon)) - 1) % nlon
        # One or two (across the 0/360 seam) contiguous longitude blocks
        blocks = [(i0, min(i0 + nx, nlon))]
        if i0 + nx > nlon:
            blocks.append((0, i0 + nx - nlon))

        def read(var, k):
            return np.concatenate([ds[var][k, a:b, :][:, j0:j1] for a, b in blocks])

        # (nc, nx, ny) in the file -> (nc, ny, nx)
        amp = np.stack([read(amp_var, k) for k in con_index]).transpose(0, 2, 1)
        phase = np.stack([read(phase_var, k) for k in con_index]).transpose(0, 2, 1)
        self.values = (amp * amp_scale * np.exp(-1j * np.radians(phase))).astype(np.complex64)
        self.lon0, self.lat0 = lon[0] + i0 * self.dlon, lat[0] + j0 * self.dlat
        self.nlat, self.nlon = self.values.shape[1], self.values.shape[2]

    def interpolate(self, lon, lat):
        """(nc, n) complex constants at lon/lat (degrees)."""
        if self.periodic:
            x = np.mod(lon - self.lon0, 360.0) / self.dlon
        else:  # positions outside the loaded block get its edge values
            half = (self.nlon - 1) * self.dlon / 2
            x = (np.mod(lon - self.lon0 - half + 180, 360.0) - 180 + half) / self.dlon
            x = np.clip(x, 0.0, self.nlon - 1.000001)
        y = np.clip((lat - self.lat0) / self.dlat, 0.0, self.nlat - 1.000001)
        i0 = np.floor(x).astype(int)
        j0 = np.floor(y).astype(int)
        wx = (x - i0).astype(np.float32)
        wy = (y - j0).astype(np.float32)
        i0 %= self.nlon
        i1 = (i0 + 1) % self.nlon
        j1 = j0 + 1
        v = self.values
        return ((1 - wy) * ((1 - wx) * v[:, j0, i0] + wx * v[:, j0, i1])
                + wy * ((1 - wx) * v[:, j1, i0] + wx * v[:, j1, i1]))


# One copy per process and set of arguments: readers rebuilt with the same
# coverage (e.g. once per simulation chunk) reuse it instead of re-reading.
_cache = {}


def _load(tpxo_dir, constituents, corners):
    """(constituents, u grid, v grid), loaded within corners and cached"""
    key = (os.path.abspath(tpxo_dir),
           None if constituents is None else tuple(constituents), tuple(corners))
    if key in _cache:
        return _cache[key]
    import netCDF4
    with netCDF4.Dataset(os.path.join(tpxo_dir, TPXO_UV_FILE)) as ds:
        ds.set_auto_mask(False)
        names = [b''.join(c).decode().strip().lower() for c in ds['con'][:]]
        if constituents is None:
            missing = [c for c in MAJOR_CONSTITUENTS if c not in names]
            if missing:
                logger.warning(f'Tidal constituents {missing} not in {TPXO_UV_FILE}, '
                               'using the others')
            constituents = tuple(c for c in MAJOR_CONSTITUENTS if c in names)
            if not constituents:
                raise ValueError(f'None of {MAJOR_CONSTITUENTS} in {TPXO_UV_FILE} ({names})')
        else:
            constituents = tuple(c.lower() for c in constituents)
            missing = [c for c in constituents if c not in names]
            if missing:
                raise ValueError(f'Constituents {missing} not in {TPXO_UV_FILE} ({names})')
        con_index = [names.index(c) for c in constituents]
        _cache[key] = (constituents,
                       HarmonicGrid(ds, 'ua', 'up', 'lon_u', 'lat_u', con_index, corners),
                       HarmonicGrid(ds, 'va', 'vp', 'lon_v', 'lat_v', con_index, corners))
    return _cache[key]


class Reader(ContinuousReader):
    """Eastward/northward TPXO barotropic tidal current (m/s).

    Args:
        tpxo_dir: folder with the TPXO OTIS netCDF velocity file (u_tpxo9.v1.nc).
        constituents: list of constituents, e.g. ['m2', 's2']. Default: those
            of the eight major ones (m2 s2 n2 k2 k1 o1 p1 q1) that the file has.
        corners: [lonmin, lonmax, latmin, latmax] to load TPXO for, e.g. the
            simulation domain. Default global (about 300 MB in memory).
    """

    def __init__(self, tpxo_dir, constituents=None, corners=(-180, 180, -90, 90)):
        try:
            import pyTMD  # noqa: F401
        except ImportError:
            raise ImportError('The pyTMD package is needed for tidal currents, '
                              'see https://pytmd.readthedocs.io')
        self.constituents, self.u, self.v = _load(tpxo_dir, constituents, corners)
        self.variables = ['x_sea_water_velocity', 'y_sea_water_velocity']
        self.proj4 = '+proj=latlong'
        self.xmin, self.xmax, self.ymin, self.ymax = corners
        self.start_time = None
        self.end_time = None
        self.time_step = None
        self.name = f'TPXO tides ({",".join(self.constituents)})'
        super().__init__()

    def get_variables(self, requested_variables, time=None, x=None, y=None, z=None):
        phasors = tidal_phasors(time, self.constituents).astype(np.complex64)
        grids = {'x_sea_water_velocity': self.u, 'y_sea_water_velocity': self.v}
        variables = {'time': time, 'x': x, 'y': y, 'z': z}
        for var in requested_variables:
            variables[var] = np.real(phasors @ grids[var].interpolate(
                np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))).astype(np.float64)
        return variables
