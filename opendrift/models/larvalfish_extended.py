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
# Copyright 2021, Trond Kristiansen, Niva
# Jan 2021 Simplified by Knut-Frode Dagestad, MET Norway, and adapted to to Kvile et al. (2018)
# 2025-2026 Extended by Kristen Thyng with configurable vertical behavior,
#           diel vertical migration (DVM), phytoplankton support, and
#           fixed-time egg hatching.

import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import Lagrangian3DArray, OceanDrift
from opendrift.models.physics_methods import solar_elevation
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED


class LarvalFishExtendedElement(Lagrangian3DArray):
    """
    Extending Lagrangian3DArray with specific properties for larval fish
    and generic biological particles (e.g. phytoplankton).
    """

    variables = Lagrangian3DArray.add_variables([
        ('stage_fraction', {'dtype': np.float32,  # to track percentage of development time completed
                            'units': '',
                            'default': 0.}),
        ('hatched', {'dtype': np.uint8,  # 0 for eggs, 1 for larvae
                     'units': '',
                     'default': 0}),
        ])


class LarvalFishExtended(OceanDrift):
    """Extended biological particle trajectory model based on OpenDrift.

    Extends the base LarvalFish model with:

    - Configurable vertical behavior: depth-keeping and diel vertical
      migration (DVM) using solar elevation for day/night detection.
    - Support for both larval fish and phytoplankton particle types.
    - Multiple egg hatching methods (temperature-dependent or fixed-time).
    - Adaptive depth-band targeting with configurable band widths.

    Vertical behavior modes
    -----------------------
    - ``none``: No active vertical movement (passive drift).
    - ``depth``: Maintain a preferred depth band around ``z_pref``.
    - ``dvm``: Diel vertical migration between ``z_day`` (daytime) and
      ``z_night`` (nighttime) depths, using ``solar_elevation()`` to
      determine day/night state.

    Particle types
    --------------
    - ``larva``: Full egg/hatching/growth lifecycle. Only hatched larvae
      exhibit active vertical behavior; eggs remain passive.
    - ``phytoplankton``: All particles have active vertical behavior
      from the start. No egg/growth stages.
    """

    ElementType = LarvalFishExtendedElement

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_height': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': 100},
        'ocean_vertical_diffusivity': {'fallback': 0.01, 'profiles': True},
        'ocean_mixed_layer_thickness': {'fallback': 50},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
    }


    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(LarvalFishExtended, self).__init__(*args, **kwargs)

        # IBM configuration options
        self._add_config({
            # Particle type selection
            'biology:particle_type':
                {'type': 'enum', 'enum': ['larva', 'phytoplankton'],
                 'default': 'larva',
                 'description': 'Type of particle to simulate. Larvae have egg/hatching/growth stages. Phytoplankton only use vertical behavior.',
                 'level': CONFIG_LEVEL_ESSENTIAL},

            # Vertical behavior system
            'biology:vertical_behavior_mode':
                {'type': 'enum', 'enum': ['none', 'depth', 'dvm'],
                 'default': 'dvm',
                 'description': 'Vertical behavior mode. none: no active movement. depth: maintain preferred depth band. dvm: diel vertical migration between day/night depths.',
                 'level': CONFIG_LEVEL_ESSENTIAL},

            'biology:w_active':
                {'type': 'float', 'default': 0.003,
                 'min': 0.0, 'max': 1.0, 'units': 'm/s',
                 'description': 'Maximum active vertical positioning speed (swimming/buoyancy regulation).',
                 'level': CONFIG_LEVEL_BASIC},

            'biology:z_pref':
                {'type': 'float', 'default': -10.0,
                 'min': -10000, 'max': 0.0, 'units': 'm',
                 'description': 'Preferred depth for depth mode (negative down from surface).',
                 'level': CONFIG_LEVEL_BASIC},

            'biology:z_day':
                {'type': 'float', 'default': -25.0,
                 'min': -10000, 'max': 0.0, 'units': 'm',
                 'description': 'Target depth during daytime for DVM mode (negative down from surface).',
                 'level': CONFIG_LEVEL_BASIC},

            'biology:z_night':
                {'type': 'float', 'default': -5.0,
                 'min': -10000, 'max': 0.0, 'units': 'm',
                 'description': 'Target depth during nighttime for DVM mode (negative down from surface).',
                 'level': CONFIG_LEVEL_BASIC},

            # Internal band expansion parameters
            'biology:dz_min':
                {'type': 'float', 'default': 1.0,
                 'min': 0.1, 'max': 100, 'units': 'm',
                 'description': 'Minimum half-width for depth bands.',
                 'level': CONFIG_LEVEL_ADVANCED},

            'biology:dz_rel':
                {'type': 'float', 'default': 0.1,
                 'min': 0.0, 'max': 1.0, 'units': 'fraction',
                 'description': 'Relative depth band expansion factor (fraction of depth).',
                 'level': CONFIG_LEVEL_ADVANCED},

            'biology:dz_max':
                {'type': 'float', 'default': 15.0,
                 'min': 0.1, 'max': 1000, 'units': 'm',
                 'description': 'Maximum half-width for depth bands.',
                 'level': CONFIG_LEVEL_ADVANCED},

            # Egg hatching options
            'egg:hatching_method':
                {'type': 'enum', 'enum': ['fixed_time'],
                 'default': 'fixed_time',
                 'description': 'Egg hatching method. fixed_time: hatch after fixed duration.',
                 'level': CONFIG_LEVEL_BASIC},

            'egg:hatch_time_days':
                {'type': 'float', 'default': 2.0,
                 'min': 0.004, 'max': 416, 'units': 'days',
                 'description': 'Fixed time to hatching when hatching_method is fixed_time.',
                 'level': CONFIG_LEVEL_BASIC},
            })

        self._set_config_default('drift:vertical_mixing', True)
        self._set_config_default('drift:vertical_mixing_at_surface', True)
        self._set_config_default('drift:vertical_advection_at_surface', True)

    # ---------------------------------------------------------------------
    # Depth band computation
    # ---------------------------------------------------------------------

    def _compute_band_half_width(self, z):
        """Compute depth band half-width: clamp(dz_min, dz_rel*|z|, dz_max)."""
        dz_min = self.get_config('biology:dz_min')
        dz_rel = self.get_config('biology:dz_rel')
        dz_max = self.get_config('biology:dz_max')

        dz = dz_rel * np.abs(z)
        dz = np.maximum(dz, dz_min)
        dz = np.minimum(dz, dz_max)
        return dz

    def _target_into_band(self, z, center, half_w):
        """Return target depth that moves z into [center-half_w, center+half_w].

        If z is already inside the band, returns z (no movement).
        """
        band_min = center - half_w
        band_max = center + half_w

        return np.where(
            (z >= band_min) & (z <= band_max),
            z,
            np.where(z < band_min, band_min, band_max),
        )

    # ---------------------------------------------------------------------
    # Vertical behavior
    # ---------------------------------------------------------------------

    def _apply_vertical_behavior(self):
        """Apply active vertical behavior based on configured mode.

        Modes
        -----
        - ``none``: No active vertical movement.
        - ``depth``: Maintain depth band around ``z_pref``.
        - ``dvm``: Diel vertical migration between ``z_day`` and ``z_night``,
          using ``solar_elevation()`` for day/night detection.

        For larvae (particle_type='larva'), only hatched larvae (not eggs)
        exhibit active vertical behavior. Eggs remain passive.
        For phytoplankton, all particles are active.
        """
        behavior = self.get_config('biology:vertical_behavior_mode')

        if behavior == 'none':
            logger.debug('_apply_vertical_behavior: mode=none, skipping')
            return

        particle_type = self.get_config('biology:particle_type')

        # Determine which particles should have active behavior
        if particle_type == 'larva':
            active_idx = np.where(self.elements.hatched == 1)[0]
            if len(active_idx) == 0:
                logger.debug('_apply_vertical_behavior: no hatched larvae, skipping')
                return
        else:
            # Phytoplankton: all particles have active behavior
            active_idx = np.arange(self.num_elements_active())

        z = self.elements.z[active_idx].copy()
        dt = self.time_step.total_seconds()
        w_active = self.get_config('biology:w_active')

        logger.debug(f'_apply_vertical_behavior: mode={behavior}, w_active={w_active}, dt={dt}s, '
                     f'active particles={len(active_idx)}, '
                     f'initial z range=[{np.min(z):.2f}, {np.max(z):.2f}]')

        if w_active <= 0.0 or dt <= 0.0:
            return

        # Determine target depth based on mode
        if behavior == 'depth':
            z_pref = self.get_config('biology:z_pref')
            half_w = self._compute_band_half_width(z_pref)
            target_z = self._target_into_band(z, z_pref, half_w)
        elif behavior == 'dvm':
            lon = np.asarray(self.elements.lon[active_idx], dtype=float)
            lat = np.asarray(self.elements.lat[active_idx], dtype=float)

            # Use OpenDrift's solar_elevation to determine day vs night
            is_day = solar_elevation(self.time, lon, lat) > 0

            z_day = self.get_config('biology:z_day')
            z_night = self.get_config('biology:z_night')
            half_w_day = self._compute_band_half_width(z_day)
            half_w_night = self._compute_band_half_width(z_night)

            target_day = self._target_into_band(z, z_day, half_w_day)
            target_night = self._target_into_band(z, z_night, half_w_night)
            target_z = np.where(is_day, target_day, target_night).astype(z.dtype)

            logger.debug(f'_apply_vertical_behavior dvm: is_day={np.sum(is_day)}/{len(is_day)}, '
                         f'z_day={z_day}, z_night={z_night}')
        else:
            logger.warning(f'Unknown vertical behavior mode: {behavior}')
            return

        # Calculate and apply displacement
        dz = target_z - z
        max_step = w_active * dt
        dz_step = np.clip(dz, -max_step, max_step)

        self.elements.z[active_idx] += dz_step

        # Clip to water column
        self.elements.z[active_idx] = np.minimum(self.elements.z[active_idx], 0.0)
        if hasattr(self.environment, 'sea_floor_depth_below_sea_level'):
            bottom = -self.environment.sea_floor_depth_below_sea_level[active_idx]
            self.elements.z[active_idx] = np.maximum(self.elements.z[active_idx], bottom)

        logger.debug(f'_apply_vertical_behavior: final z range=[{np.min(self.elements.z[active_idx]):.2f}, '
                     f'{np.max(self.elements.z[active_idx]):.2f}]')

    def update_fish_larvae(self):

        # Hatching of eggs
        eggs = np.where(self.elements.hatched==0)[0]
        if len(eggs) > 0:
            hatching_method = self.get_config('egg:hatching_method')

            if hatching_method == 'fixed_time':
                # Fixed-time hatching
                hatch_time_days = self.get_config('egg:hatch_time_days')
                days_in_timestep = self.time_step.total_seconds() / 86400
                self.elements.stage_fraction[eggs] += days_in_timestep / hatch_time_days
                hatching = np.where(self.elements.stage_fraction[eggs]>=1)[0]

            else:
                logger.warning(f'Unknown hatching method: {hatching_method}')
                hatching = []

            if len(hatching) > 0:
                logger.debug('Hatching %s eggs' % len(hatching))
                self.elements.hatched[eggs[hatching]] = 1

        larvae = np.where(self.elements.hatched==1)[0]
        if len(larvae) == 0:
            logger.debug('%s eggs, with maximum stage_fraction of %s (1 gives hatching)'
                         % (len(eggs), self.elements.stage_fraction[eggs].max()))
            return

    # ---------------------------------------------------------------------
    # Main update loop
    # ---------------------------------------------------------------------

    def update(self):
        """Update particle state: biology, advection, mixing, vertical behavior."""
        particle_type = self.get_config('biology:particle_type')

        # Fish-specific processes (eggs, larvae, growth)
        if particle_type == 'larva':
            self.update_fish_larvae()

        # Physical advection
        self.advect_ocean_current()

        # Stokes drift
        self.stokes_drift()

        # Vertical mixing
        self.vertical_mixing()

        # Active vertical behavior (depth-keeping or DVM)
        self._apply_vertical_behavior()
