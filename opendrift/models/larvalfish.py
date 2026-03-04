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

import datetime
import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import Lagrangian3DArray, OceanDrift
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED


class LarvalFishElement(Lagrangian3DArray):
    """
    Extending Lagrangian3DArray with specific properties for larval and juvenile stages of fish
    """

    variables = Lagrangian3DArray.add_variables([
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0014}),  # for NEA Cod
        ('neutral_buoyancy_salinity', {'dtype': np.float32,
                                       'units': 'PSU',
                                       'default': 31.25}),  # for NEA Cod
        ('stage_fraction', {'dtype': np.float32,  # to track percentage of development time completed
                            'units': '',
                            'default': 0.}),
        ('hatched', {'dtype': np.uint8,  # 0 for eggs, 1 for larvae
                     'units': '',
                     'default': 0}),
        ('length', {'dtype': np.float32,
                    'units': 'mm',
                    'default': 0}),
        ('weight', {'dtype': np.float32,
                    'units': 'mg',
                    'default': 0.08}),
        ('survival', {'dtype': np.float32,  # Not yet used
                      'units': '',
                      'default': 1.})])


class LarvalFish(OceanDrift):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

    """

    ElementType = LarvalFishElement

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
        super(LarvalFish, self).__init__(*args, **kwargs)

        # IBM configuration options
        self._add_config({
            'IBM:fraction_of_timestep_swimming':
                {'type': 'float', 'default': 0.15,
                 'min': 0.0, 'max': 1.0, 'units': 'fraction',
                 'description': 'Fraction of timestep swimming',
                 'level': CONFIG_LEVEL_ADVANCED},
            
            # Particle type selection
            'biology:particle_type':
                {'type': 'enum', 'enum': ['larva', 'phytoplankton'],
                 'default': 'larva',
                 'description': 'Type of particle to simulate. Larvae have egg/hatching/growth stages. Phytoplankton only use vertical behavior.',
                 'level': CONFIG_LEVEL_ESSENTIAL},
            
            # Unified vertical behavior system
            'biology:vertical_behavior_mode':
                {'type': 'enum', 'enum': ['none', 'depth', 'dvm', 'legacy'],
                 'default': 'legacy',
                 'description': 'Vertical behavior mode. none: no active movement. depth: maintain preferred depth band. dvm: diel vertical migration between day/night depths. legacy: original larvalfish time-based swimming.',
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
            
            # Internal band expansion parameters (global)
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
                {'type': 'enum', 'enum': ['temperature', 'fixed_time'],
                 'default': 'temperature',
                 'description': 'Egg hatching method. temperature: use ambient temperature (Ellertsen et al. 1988). fixed_time: hatch after fixed duration.',
                 'level': CONFIG_LEVEL_BASIC},
            
            'egg:hatch_time_hours':
                {'type': 'float', 'default': 48.0,
                 'min': 0.1, 'max': 10000, 'units': 'hours',
                 'description': 'Fixed time to hatching when hatching_method is fixed_time.',
                 'level': CONFIG_LEVEL_BASIC},
            })

        self._set_config_default('drift:vertical_mixing', True)
        self._set_config_default('drift:vertical_mixing_at_surface', True)
        self._set_config_default('drift:vertical_advection_at_surface', True)

    # ---------------------------------------------------------------------
    # Time helpers for DVM: UTC -> local solar hour
    # ---------------------------------------------------------------------

    def _current_utc_hour(self):
        """Current model time in UTC hours [0, 24)."""
        t = self.time  # datetime in UTC
        return t.hour + t.minute / 60.0 + t.second / 3600.0

    def _local_solar_hour(self, lon_deg):
        """
        Convert UTC model time to local solar time for given longitude(s).
        
        Parameters
        ----------
        lon_deg : float or array
            Longitude in degrees East (scalar or array).
        
        Returns
        -------
        float or array
            Local solar hours in [0, 24).
        """
        utc_hour = self._current_utc_hour()
        # 15° longitude = 1 hour; positive lon east of Greenwich.
        local = utc_hour + lon_deg / 15.0
        # Wrap into [0, 24)
        return np.mod(local, 24.0)

    def _approx_sunrise_sunset(self, lat_deg):
        """
        Fast approximation of sunrise/sunset for current simulation date.
        
        Uses day-of-year from self.time (UTC) to estimate sunrise/sunset times
        based on latitude. This is a simplified calculation suitable for 
        particle tracking applications.
        
        Parameters
        ----------
        lat_deg : float or array
            Latitude in degrees (scalar or array).
        
        Returns
        -------
        sunrise : float or array
            Sunrise time in local solar hours [0, 24).
        sunset : float or array
            Sunset time in local solar hours [0, 24).
        
        References
        ----------
        Simplified from various solar position algorithms. Accuracy is 
        sufficient for DVM behavior modeling.
        """
        # Day of year from the simulation clock (UTC)
        N = self.time.timetuple().tm_yday

        # Solar declination (degrees) - approximation
        decl_deg = 23.44 * np.sin(2 * np.pi * (N - 80) / 365)
        decl_rad = np.radians(decl_deg)
        lat_rad = np.radians(lat_deg)

        # Hour angle at sunrise/sunset
        cosH = -np.tan(lat_rad) * np.tan(decl_rad)
        cosH = np.clip(cosH, -1.0, 1.0)   # handle polar extremes
        H = np.arccos(cosH)               # radians

        # Convert to hours: 15° per hour
        H_deg = np.degrees(H)
        daylight_hours = 2.0 * H_deg / 15.0

        sunrise = 12.0 - daylight_hours / 2.0
        sunset = 12.0 + daylight_hours / 2.0

        return sunrise, sunset

    # ---------------------------------------------------------------------
    # Depth band computation with internal expansion
    # ---------------------------------------------------------------------

    def _compute_band_half_width(self, z):
        """
        Compute depth band half-width using adaptive expansion rule.
        
        Parameters
        ----------
        z : float or array
            Target depth (negative down from surface, in meters).
        
        Returns
        -------
        float or array
            Half-width of depth band in meters.
        
        Notes
        -----
        Uses rule: dz = clamp(dz_min, dz_rel * |z|, dz_max)
        where dz_min, dz_rel, dz_max are config parameters.
        """
        dz_min = self.get_config('biology:dz_min')
        dz_rel = self.get_config('biology:dz_rel')
        dz_max = self.get_config('biology:dz_max')
        
        # Calculate relative band width
        dz = dz_rel * np.abs(z)
        
        # Apply min/max bounds
        dz = np.maximum(dz, dz_min)
        dz = np.minimum(dz, dz_max)
        
        return dz

    def _target_into_band(self, z, center, half_w):
        """
        Return target depth that moves z into [center-half_w, center+half_w].
        
        Parameters
        ----------
        z : array
            Current particle depths.
        center : float
            Band center depth (negative down).
        half_w : float or array
            Band half-width.
        
        Returns
        -------
        array
            Target depths. If z is inside band, returns z (no movement).
            If outside band, returns nearest band edge.
        """
        band_min = center - half_w
        band_max = center + half_w

        return np.where(
            (z >= band_min) & (z <= band_max),
            z,
            np.where(z < band_min, band_min, band_max),
        )

    def _clip_to_water_column(self, z):
        """
        Clip depths to physical water column (surface to seafloor).
        
        Parameters
        ----------
        z : array
            Particle depths (negative down from surface).
        
        Returns
        -------
        array
            Clipped depths within [seafloor, 0].
        """
        # Clip to surface (0 m)
        z_clipped = np.minimum(z, 0.0)
        
        # Clip to seafloor if available
        if hasattr(self.environment, 'sea_floor_depth_below_sea_level'):
            bottom = -self.environment.sea_floor_depth_below_sea_level
            z_clipped = np.maximum(z_clipped, bottom)
        
        return z_clipped

    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        """Calculate terminal velocity for Pelagic Egg

        according to
        S. Sundby (1983): A one-dimensional model for the vertical
        distribution of pelagic fish eggs in the mixed layer
        Deep Sea Research (30) pp. 645-661

        Method copied from ibm.f90 module of LADIM:
        Vikebo, F., S. Sundby, B. Aadlandsvik and O. Otteraa (2007),
        Fish. Oceanogr. (16) pp. 216-228
        """
        g = 9.81  # ms-2

        # Pelagic Egg properties that determine buoyancy
        eggsize = self.elements.diameter  # 0.0014 for NEA Cod
        eggsalinity = self.elements.neutral_buoyancy_salinity
        # 31.25 for NEA Cod

        # prepare interpolation of temp, salt
        if not (Tprofiles is None and Sprofiles is None):
            if z_index is None:
                z_i = range(Tprofiles.shape[0])  # evtl. move out of loop
                # evtl. move out of loop
                z_index = interp1d(-self.environment_profiles['z'],
                                   z_i, bounds_error=False)
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.uint8), 0)
            lower = np.minimum(upper + 1, Tprofiles.shape[0] - 1)
            weight_upper = 1 - (zi - upper)

        # do interpolation of temp, salt if profiles were passed into
        # this function, if not, use reader by calling self.environment
        if Tprofiles is None:
            T0 = self.environment.sea_water_temperature
        else:
            T0 = Tprofiles[upper, range(Tprofiles.shape[1])] * \
                 weight_upper + \
                 Tprofiles[lower, range(Tprofiles.shape[1])] * \
                 (1 - weight_upper)
        if Sprofiles is None:
            S0 = self.environment.sea_water_salinity
        else:
            S0 = Sprofiles[upper, range(Sprofiles.shape[1])] * \
                 weight_upper + \
                 Sprofiles[lower, range(Sprofiles.shape[1])] * \
                 (1 - weight_upper)

        # The density difference between a pelagic egg and the ambient water
        # is regulated by their salinity difference through the
        # equation of state for sea water.
        # The Egg has the same temperature as the ambient water and its
        # salinity is regulated by osmosis through the egg shell.
        DENSw = self.sea_water_density(T=T0, S=S0)
        DENSegg = self.sea_water_density(T=T0, S=eggsalinity)
        dr = DENSw - DENSegg  # density difference

        # water viscosity
        my_w = 0.001 * (1.7915 - 0.0538 * T0 + 0.007 * (T0 ** (2.0)) - 0.0023 * S0)
        # ~0.0014 kg m-1 s-1

        # terminal velocity for low Reynolds numbers
        W = (1.0 / my_w) * (1.0 / 18.0) * g * eggsize ** 2 * dr

        # check if we are in a Reynolds regime where Re > 0.5
        highRe = np.where(W * 1000 * eggsize / my_w > 0.5)

        # Use empirical equations for terminal velocity in
        # high Reynolds numbers.
        # Empirical equations have length units in cm!
        my_w = 0.01854 * np.exp(-0.02783 * T0)  # in cm2/s
        d0 = (eggsize * 100) - 0.4 * \
             (9.0 * my_w ** 2 / (100 * g) * DENSw / dr) ** (1.0 / 3.0)  # cm
        W2 = 19.0 * d0 * (0.001 * dr) ** (2.0 / 3.0) * (my_w * 0.001 * DENSw) ** (-1.0 / 3.0)
        # cm/s
        W2 = W2 / 100.  # back to m/s

        W[highRe] = W2[highRe]
        self.elements.terminal_velocity = W

    # ---------------------------------------------------------------------
    # Vertical behavior targeting
    # ---------------------------------------------------------------------

    def _target_depth_mode(self, z):
        """
        Target depth for 'depth' behavior mode: maintain preferred depth band.
        
        Parameters
        ----------
        z : array
            Current particle depths.
        
        Returns
        -------
        array
            Target depths. Inside band: no change. Outside: swim to band edge.
        """
        z_pref = self.get_config('biology:z_pref')
        half_w = self._compute_band_half_width(z_pref)
        
        target = self._target_into_band(z, z_pref, half_w)
        
        logger.debug(f'_target_depth_mode: z_pref={z_pref}, half_w={half_w:.2f}, '
                     f'z range=[{np.min(z):.2f}, {np.max(z):.2f}], '
                     f'target range=[{np.min(target):.2f}, {np.max(target):.2f}]')
        
        return target

    def _target_depth_dvm(self, z):
        """
        Target depth for 'dvm' (diel vertical migration) behavior mode.
        
        Switches between day and night depth bands based on local solar time
        and sunrise/sunset calculation.
        
        Parameters
        ----------
        z : array
            Current particle depths.
        
        Returns
        -------
        array
            Target depths based on day/night state for each particle.
        """
        lat = np.asarray(self.elements.lat, dtype=float)
        lon = np.asarray(self.elements.lon, dtype=float)

        # Local solar hour for each particle
        local_hour = self._local_solar_hour(lon)

        # Sunrise/sunset (local solar time) based on lat + date
        sunrise, sunset = self._approx_sunrise_sunset(lat)

        # Daytime mask
        is_day = (local_hour >= sunrise) & (local_hour < sunset)

        z_day = self.get_config('biology:z_day')
        z_night = self.get_config('biology:z_night')
        
        # Debug logging
        utc_hour = self._current_utc_hour()
        logger.debug(f'_target_depth_dvm: UTC hour={utc_hour:.1f}, '
                     f'local_hour={np.mean(local_hour):.1f}, '
                     f'sunrise={np.mean(sunrise):.1f}, sunset={np.mean(sunset):.1f}, '
                     f'is_day={np.sum(is_day)}/{len(is_day)} particles, '
                     f'z_day={z_day}, z_night={z_night}')
        
        # Compute adaptive band widths for each target depth
        half_w_day = self._compute_band_half_width(z_day)
        half_w_night = self._compute_band_half_width(z_night)

        # Targets for day and night bands
        target_day = self._target_into_band(z, z_day, half_w_day)
        target_night = self._target_into_band(z, z_night, half_w_night)

        # Choose per-particle target based on day/night
        target = np.where(is_day, target_day, target_night).astype(z.dtype)
        return target

    def _apply_vertical_behavior(self):
        """
        Apply active vertical behavior based on configured mode.
        
        Implements unified vertical behavior system for both larvae and 
        phytoplankton. Particles swim/regulate buoyancy toward target depths
        with speed limited by w_active parameter.
        
        For larvae (particle_type='larva'), only hatched larvae (not eggs)
        exhibit active vertical behavior. Eggs remain passive.
        
        Notes
        -----
        Modes:
        - 'none': No active vertical movement
        - 'depth': Maintain depth band around z_pref
        - 'dvm': Diel vertical migration between z_day and z_night
        - 'legacy': Use original larvae_vertical_migration method
        """
        behavior = self.get_config('biology:vertical_behavior_mode')
        
        if behavior == 'none':
            logger.debug('_apply_vertical_behavior: mode=none, skipping')
            return
        
        if behavior == 'legacy':
            logger.debug('_apply_vertical_behavior: mode=legacy, skipping (called separately)')
            return
        
        particle_type = self.get_config('biology:particle_type')
        
        # Determine which particles should have active behavior
        if particle_type == 'larva':
            # Only hatched larvae (not eggs) have active vertical behavior
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

        # Get lat/lon for active particles only
        lat = np.asarray(self.elements.lat[active_idx], dtype=float)
        lon = np.asarray(self.elements.lon[active_idx], dtype=float)

        # Determine target depth based on mode
        if behavior == 'depth':
            z_pref = self.get_config('biology:z_pref')
            half_w = self._compute_band_half_width(z_pref)
            target_z = self._target_into_band(z, z_pref, half_w)
        elif behavior == 'dvm':
            local_hour = self._local_solar_hour(lon)
            sunrise, sunset = self._approx_sunrise_sunset(lat)
            is_day = (local_hour >= sunrise) & (local_hour < sunset)
            
            z_day = self.get_config('biology:z_day')
            z_night = self.get_config('biology:z_night')
            half_w_day = self._compute_band_half_width(z_day)
            half_w_night = self._compute_band_half_width(z_night)
            
            target_day = self._target_into_band(z, z_day, half_w_day)
            target_night = self._target_into_band(z, z_night, half_w_night)
            target_z = np.where(is_day, target_day, target_night).astype(z.dtype)
            
            logger.debug(f'_target_depth_dvm: is_day={np.sum(is_day)}/{len(is_day)}, '
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

    def fish_growth(self, weight, temperature):
        # Weight in milligrams, temperature in celcius
        # Daily growth rate in percent according to
        # Folkvord, A. 2005. "Comparison of Size-at-Age of Larval Atlantic Cod (Gadus Morhua)
        # from Different Populations Based on Size- and Temperature-Dependent Growth
        # Models." Canadian Journal of Fisheries and Aquatic Sciences.
        # Journal Canadien Des Sciences Halieutiques # et Aquatiques 62(5): 1037-52.
        GR = 1.08 + 1.79 * temperature - 0.074 * temperature * np.log(weight) \
             - 0.0965 * temperature * np.log(weight) ** 2 \
             + 0.0112 * temperature * np.log(weight) ** 3

        # Growth rate(g) converted to milligram weight (gr_mg) per timestep:
        g = (np.log(GR / 100. + 1)) * self.time_step.total_seconds()/86400
        return weight * (np.exp(g) - 1.)

    def update_fish_larvae(self):

        # Hatching of eggs
        eggs = np.where(self.elements.hatched==0)[0]
        if len(eggs) > 0:
            hatching_method = self.get_config('egg:hatching_method')
            
            if hatching_method == 'temperature':
                # Original temperature-dependent hatching (Ellertsen et al. 1988)
                amb_duration = np.exp(3.65 - 0.145*self.environment.sea_water_temperature[eggs]) # Total egg development time (days) according to ambient temperature
                days_in_timestep = self.time_step.total_seconds()/(60*60*24)  # The fraction of a day completed in one time step
                amb_fraction = days_in_timestep/amb_duration # Fraction of development time completed during present time step
                self.elements.stage_fraction[eggs] += amb_fraction # Add fraction completed during present timestep to cumulative fraction completed
                hatching = np.where(self.elements.stage_fraction[eggs]>=1)[0]
                
            elif hatching_method == 'fixed_time':
                # Fixed-time hatching
                hatch_time_hours = self.get_config('egg:hatch_time_hours')
                hatch_time_seconds = hatch_time_hours * 3600
                hours_in_timestep = self.time_step.total_seconds() / 3600
                # Use stage_fraction to track time: fraction = (time elapsed) / (hatch_time)
                self.elements.stage_fraction[eggs] += hours_in_timestep / hatch_time_hours
                hatching = np.where(self.elements.stage_fraction[eggs]>=1)[0]
                
            else:
                logger.warning(f'Unknown hatching method: {hatching_method}')
                hatching = []
            
            if len(hatching) > 0:
                logger.debug('Hatching %s eggs' % len(hatching))
                self.elements.hatched[eggs[hatching]] = 1 # Eggs with total development time completed are hatched (1)

        larvae = np.where(self.elements.hatched==1)[0]
        if len(larvae) == 0:
            logger.debug('%s eggs, with maximum stage_fraction of %s (1 gives hatching)'
                         % (len(eggs), self.elements.stage_fraction[eggs].max()))
            return

        # Increasing weight of larvae
        avg_weight_before = self.elements.weight[larvae].mean()
        growth = self.fish_growth(self.elements.weight[larvae],
                                  self.environment.sea_water_temperature[larvae])
        self.elements.weight[larvae] += growth
        avg_weight_after = self.elements.weight[larvae].mean()
        logger.debug('Growing %s larve from average size %s to %s' %
              (len(larvae), avg_weight_before, avg_weight_after))

        # Increasing length of larvae, according to Folkvord (2005)
        w = self.elements.weight[larvae]
        self.elements.length[larvae] = np.exp(2.296 + 0.277 * np.log(w) - 0.005128 *np.log10(w)**2)

    def larvae_vertical_migration(self):
        """
        Original legacy vertical migration behavior.
        
        Simple time-of-day based swimming: down before noon, up after noon.
        Only used when biology:vertical_behavior_mode is 'legacy'.
        
        References
        ----------
        Based on Peck et al. 2006 swimming speed formula.
        """
        behavior = self.get_config('biology:vertical_behavior_mode')
        
        if behavior != 'legacy':
            return  # New modes handled by _apply_vertical_behavior
        
        larvae = np.where(self.elements.hatched==1)[0]
        if len(larvae) == 0:
            return

        # Vertical migration of Larvae
        # Swim function from Peck et al. 2006
        L = self.elements.length[larvae]
        swim_speed = (0.261*(L**(1.552*L**(-0.08))) - 5.289/L) / 1000
        f = self.get_config('IBM:fraction_of_timestep_swimming')
        max_migration_per_timestep = f*swim_speed*self.time_step.total_seconds()

        # Using here UTC hours. Should be changed to local solar time,
        # although a phase shift of some hours should not make much difference
        if self.time.hour < 12:
            direction = -1  # Swimming down when light is increasing
        else:
            direction = 1  # Swimming up when light is decreasing

        self.elements.z[larvae] = np.minimum(0, self.elements.z[larvae] + direction*max_migration_per_timestep)

    def update(self):
        """
        Main update loop for particle advection and behavior.
        
        Conditionally applies fish-specific processes (egg hatching, growth)
        based on particle_type configuration. Always applies vertical behavior
        and physical transport.
        """
        particle_type = self.get_config('biology:particle_type')
        behavior_mode = self.get_config('biology:vertical_behavior_mode')
        
        # Fish-specific processes (eggs, larvae, growth)
        if particle_type == 'larva':
            self.update_fish_larvae()
        
        # Physical advection
        self.advect_ocean_current()

        # Stokes drift
        self.stokes_drift()

        # Buoyancy (for eggs)
        if behavior_mode == 'legacy' and particle_type == 'larva':
            self.update_terminal_velocity()
        
        # Vertical mixing
        self.vertical_mixing()
        
        # Vertical behavior
        if behavior_mode == 'legacy' and particle_type == 'larva':
            # Legacy mode: only for larvae, uses original method
            self.larvae_vertical_migration()
        elif behavior_mode in ['depth', 'dvm']:
            # New unified behavior: works for both larvae and phytoplankton
            self._apply_vertical_behavior()
        # else: behavior_mode == 'none', no active vertical movement
