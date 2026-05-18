#!/usr/bin/env python
"""
AMS → BGO contrail drift
========================

Contrail segments are seeded continuously along a flight from Amsterdam
Schiphol (AMS) to Bergen Airport Flesland (BGO) at cruise altitude
(10 500 m). The segments then drift eastward on the upper-tropospheric
jet stream and slowly sublimate. A ✈ icon follows the flight route.
"""

import numpy as np
from datetime import datetime, timedelta
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as mplanim
from opendrift.models.contrail import ContrailDrift

# ── Flight parameters ─────────────────────────────────────────────────────────
LON_AMS, LAT_AMS = 4.76,  52.31   # Amsterdam Schiphol (AMS)
LON_BGN, LAT_BGN = 5.22,  60.29   # Bergen Airport Flesland (BGO)
CRUISE_ALT_M     = 10_500          # cruise altitude [m]
T_DEPART         = datetime(2024, 6, 15, 10, 0)   # UTC
FLIGHT_DUR       = timedelta(minutes=90)

# ── Simulation parameters ─────────────────────────────────────────────────────
SIM_DURATION = timedelta(hours=6)
TIME_STEP    = 300     # calculation step [s]
OUTPUT_STEP  = 600     # output / animation step [s]
N_SEEDS      = 18      # seeding events along route (≈ every 5 minutes)
N_PER_SEED   = 4       # elements per seeding event

# ── Build model with constant upper-tropospheric forcing ──────────────────────
c = ContrailDrift(loglevel=30)

# Constant environment – replace with a real NWP reader for a realistic run
c.set_config('environment:constant:x_wind',             22)     # jet stream [m s-1]
c.set_config('environment:constant:y_wind',              2)     # slight northward
c.set_config('environment:constant:air_temperature',   220)     # [K]
c.set_config('environment:constant:specific_humidity', 8.5e-5)  # near ice-saturation [kg kg-1]
c.set_config('environment:constant:horizontal_diffusivity', 80) # [m2 s-1]
c.set_config('contrail:sublimation_timescale',        5400)     # 1.5-h lifetime [s]

# ── Seed contrail segments continuously along the flight path ─────────────────
fracs      = np.linspace(0, 1, N_SEEDS)
lons_seed  = LON_AMS + fracs * (LON_BGN - LON_AMS)
lats_seed  = LAT_AMS + fracs * (LAT_BGN - LAT_AMS)
times_seed = [T_DEPART + timedelta(seconds=float(f) * FLIGHT_DUR.total_seconds())
              for f in fracs]

for lon, lat, t in zip(lons_seed, lats_seed, times_seed):
    c.seed_elements(lon=lon, lat=lat, z=CRUISE_ALT_M,
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
                       f'Contrail drift (constant 22 m s\u207b\u00b9 westerly)')

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
       title='AMS → BGO contrail drift (6 h)')
