#!/usr/bin/env python
"""
AMS → BGO contrail drift with NOAA GFS
=======================================

Contrail segments are seeded continuously along a flight from Amsterdam
Schiphol (AMS) to Bergen Airport Flesland (BGO) at **5 km altitude**.
Wind and temperature are taken from the NOAA GFS forecast (UCAR Thredds
OPeNDAP). Because GFS stores data on pressure levels, the isobaric level
nearest to 5 km (~550 hPa) is pre-selected with xarray before handing
the 2-D field to the OpenDrift reader. A ✈ icon follows the flight route.
"""

import numpy as np
from datetime import datetime, timedelta
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.animation as mplanim
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.contrail import ContrailDrift

# ── Flight parameters ─────────────────────────────────────────────────────────
LON_AMS, LAT_AMS = 4.76,  52.31   # Amsterdam Schiphol (AMS)
LON_BGN, LAT_BGN = 5.22,  60.29   # Bergen Airport Flesland (BGO)
ALTITUDE_M       = 5_000           # 5 km ≈ 550 hPa
# Use current UTC time so the GFS Best dataset is always valid
T_DEPART  = datetime.utcnow().replace(minute=0, second=0, microsecond=0)
FLIGHT_DUR = timedelta(minutes=90)

# ── Load NOAA GFS at 5 km from UCAR Thredds OPeNDAP ──────────────────────────
# The GFS "Best" dataset always reflects the most recent available forecast.
# Wind and temperature are stored on isobaric (pressure) levels in Pa.
# We pre-select the level closest to 5 km so the reader sees a simple
# 2-D (time × lat × lon) field without a vertical dimension.
GFS_URL = ('https://thredds.ucar.edu/thredds/dodsC/'
           'grib/NCEP/GFS/Global_0p25deg/Best')

# ISA pressure at 5 km: p = 101325 * (1 - 2.2557e-5 * z)^5.2559 ≈ 54 050 Pa
_isa_pa = 101325.0 * (1.0 - 2.2557e-5 * ALTITUDE_M) ** 5.2559

print(f'Opening GFS dataset: {GFS_URL}')
ds_gfs = xr.open_dataset(GFS_URL, engine='netcdf4')


def _select_level(ds, var_name, target_pa):
    """Return the DataArray for *var_name* at the isobaric level nearest to
    *target_pa* (Pascals), handling both Pa and hPa coordinate units."""
    var = ds[var_name]
    # Identify the pressure dimension (not time / lat / lon)
    iso_dim = next(
        d for d in var.dims
        if d not in {'time', 'time1', 'lat', 'lon', 'latitude', 'longitude'})
    iso_vals = ds[iso_dim].values.astype(float)
    # Detect Pa vs hPa by magnitude
    target = target_pa if iso_vals.max() > 2000.0 else target_pa / 100.0
    idx = int(np.argmin(np.abs(iso_vals - target)))
    print(f'  {var_name}: selecting {iso_dim}={iso_vals[idx]:.0f} '
          f'({"Pa" if iso_vals.max() > 2000 else "hPa"})')
    return var.isel({iso_dim: idx})


# Assemble a 2-D (time, lat, lon) xarray Dataset at the chosen level.
# Variable names that carry CF standard_name attributes are recognised and
# mapped automatically by reader_netCDF_CF_generic:
#   eastward_wind  → x_wind
#   northward_wind → y_wind
#   air_temperature → air_temperature
#   specific_humidity → specific_humidity
gfs_vars = {
    'u-component_of_wind_isobaric': 'u-component_of_wind_isobaric',
    'v-component_of_wind_isobaric': 'v-component_of_wind_isobaric',
    'Temperature_isobaric':         'Temperature_isobaric',
}
# specific_humidity is not always available; try to include it
try:
    gfs_vars['Specific_humidity_isobaric'] = 'Specific_humidity_isobaric'
    _ = ds_gfs['Specific_humidity_isobaric']
except KeyError:
    print('  Specific_humidity_isobaric not found in dataset – using fallback')
    gfs_vars.pop('Specific_humidity_isobaric', None)

ds_level = xr.Dataset(
    {name: _select_level(ds_gfs, name, _isa_pa) for name in gfs_vars})

reader_gfs = reader_netCDF_CF_generic.Reader(ds_level)
print(reader_gfs)

# ── Simulation parameters ─────────────────────────────────────────────────────
SIM_DURATION = timedelta(hours=4)
TIME_STEP    = 300     # calculation step [s]
OUTPUT_STEP  = 600     # output / animation step [s]
N_SEEDS      = 18      # seeding events along route (≈ every 5 minutes)
N_PER_SEED   = 4       # elements per seeding event

# ── Build contrail model ───────────────────────────────────────────────────────
c = ContrailDrift(loglevel=20)
c.add_reader(reader_gfs)

# Fallback values for variables not provided by this GFS level
c.set_config('environment:fallback:specific_humidity',      2e-4)   # kg kg-1
c.set_config('environment:fallback:horizontal_diffusivity',   50)   # m2 s-1
# At 5 km the air is typically subsaturated → contrails are short-lived
c.set_config('contrail:sublimation_timescale', 1800)   # 30-min timescale [s]

# ── Seed contrail segments continuously along the flight path ─────────────────
fracs      = np.linspace(0, 1, N_SEEDS)
lons_seed  = LON_AMS + fracs * (LON_BGN - LON_AMS)
lats_seed  = LAT_AMS + fracs * (LAT_BGN - LAT_AMS)
times_seed = [T_DEPART + timedelta(seconds=float(f) * FLIGHT_DUR.total_seconds())
              for f in fracs]

for lon, lat, t in zip(lons_seed, lats_seed, times_seed):
    c.seed_elements(lon=lon, lat=lat, z=ALTITUDE_M,
                    radius=2_000, number=N_PER_SEED, time=t)

# ── Run the simulation ────────────────────────────────────────────────────────
c.run(duration=SIM_DURATION, time_step=TIME_STEP, time_step_output=OUTPUT_STEP)

# ── Derive plane position at each output time step ────────────────────────────
out_times = c.result.time.values          # numpy datetime64 array
t0_ns     = np.datetime64(T_DEPART)
tend_ns   = np.datetime64(T_DEPART + FLIGHT_DUR)
flight_s  = FLIGHT_DUR.total_seconds()

plane_lon = np.full(len(out_times), np.nan)
plane_lat = np.full(len(out_times), np.nan)
for k, t in enumerate(out_times):
    elapsed_s = float((t - t0_ns) / np.timedelta64(1, 's'))
    if 0 <= elapsed_s <= flight_s:
        frac = elapsed_s / flight_s
        plane_lon[k] = LON_AMS + frac * (LON_BGN - LON_AMS)
        plane_lat[k] = LAT_AMS + frac * (LAT_BGN - LAT_AMS)

# ── Build the map with OpenDrift's helper ─────────────────────────────────────
corners = [LON_AMS - 0.5, LON_BGN + 8.0,     # extra room eastward for drift
           LAT_AMS - 0.8, LAT_BGN + 0.8]

fig, ax, _crs_plot, _x, _y, _i0, _i1 = c.set_up_map(
    corners=corners,
    ocean_color='aliceblue',
    land_color='wheat',
    lscale='i')

crs_ll = c.crs_lonlat   # PlateCarree – used for lon/lat coordinate artists

# Dashed great-circle guide line
ax.plot([LON_AMS, LON_BGN], [LAT_AMS, LAT_BGN],
        '--', color='gray', lw=0.8, alpha=0.5, transform=crs_ll, zorder=3)

# Airport markers and labels
for lon, lat, lbl, xoff in [(LON_AMS, LAT_AMS, 'AMS', 0.1),
                              (LON_BGN, LAT_BGN, 'BGO', 0.1)]:
    ax.plot(lon, lat, 'o', ms=5, color='steelblue', transform=crs_ll, zorder=5)
    ax.text(lon + xoff, lat - 0.2, lbl, transform=crs_ll,
            fontsize=8, color='steelblue', fontweight='bold', zorder=5)

# ── Animated artists ──────────────────────────────────────────────────────────
# Active contrail elements (cyan) and sublimated remnants (light gray)
scat_active = ax.scatter([], [], s=14, c='cyan', alpha=0.8, zorder=6,
                          transform=crs_ll, label='Active contrail')
scat_faded  = ax.scatter([], [], s=7,  c='lightgray', alpha=0.3, zorder=5,
                          transform=crs_ll, label='Sublimated')

# Plane icon – the ✈ character (U+2708) in DejaVu Sans
plane_icon = ax.text(LON_AMS, LAT_AMS, '\u2708',
                     transform=crs_ll,
                     fontsize=22, ha='center', va='center',
                     rotation=90,            # nose pointing north
                     color='navy', zorder=8,
                     visible=False)

title_obj = ax.set_title('')
ax.legend(loc='lower right', fontsize=8, framealpha=0.85)

# ── Pre-extract trajectory arrays ─────────────────────────────────────────────
lons_t = c.result.lon.values     # (n_traj, n_time)
lats_t = c.result.lat.values
n_times = lons_t.shape[1]

def animate(i):
    t_str = str(np.datetime_as_string(out_times[i], unit='m')).replace('T', ' ')
    title_obj.set_text(f'AMS \u2192 BGO  \u2502  {t_str} UTC\n'
                       f'Contrail drift at 5 km (NOAA GFS wind)')

    # Active particles: not NaN at this time step
    active = ~np.isnan(lons_t[:, i])
    faded  = np.zeros(lons_t.shape[0], dtype=bool)

    # Mark as faded: was seeded (not NaN at some earlier step) but now NaN
    if i > 0:
        was_active = ~np.isnan(lons_t[:, max(i - 1, 0)])
        faded = was_active & ~active

    scat_active.set_offsets(
        np.c_[lons_t[active, i], lats_t[active, i]]
        if active.any() else np.empty((0, 2)))

    scat_faded.set_offsets(
        np.c_[lons_t[faded, max(i - 1, 0)], lats_t[faded, max(i - 1, 0)]]
        if faded.any() else np.empty((0, 2)))

    # Plane icon: visible only while airborne
    if not np.isnan(plane_lon[i]):
        plane_icon.set_position((plane_lon[i], plane_lat[i]))
        plane_icon.set_visible(True)
    else:
        plane_icon.set_visible(False)

    return scat_active, scat_faded, plane_icon, title_obj


ani = mplanim.FuncAnimation(
    fig, animate, frames=n_times, interval=150, blit=False)

#%%
# Save animation – change filename to None to display interactively instead.
ani.save('example_contrail_ams_bgo.gif', writer='pillow', fps=6,
         savefig_kwargs={'facecolor': fig.get_facecolor()})

#%%
# .. image:: /gallery/animations/example_contrail_ams_bgo_0.gif

#%%
# Static plot of all trajectories
c.plot(fast=True, ocean_color='aliceblue', land_color='wheat',
       corners=corners, show_elements=True,
       title='AMS → BGO contrail drift at 5 km – NOAA GFS')
