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
# Copyright 2024, OpenDrift Developers

"""Lagrangian model for aircraft contrail transport and microphysics.

Each particle represents a small segment of a contrail plume in the upper
troposphere.  The model resolves:

* **Schmidt-Appleman Criterion (SAC)**: elements are deactivated immediately
  when ambient temperature exceeds the critical temperature T_SAC, i.e. when
  conditions do not support contrail formation.
* **Horizontal advection** by upper-tropospheric winds.
* **Horizontal diffusion** due to unresolved wind shear and turbulence.
* **Gravitational sedimentation** of ice crystals.
* **Sublimation** when ambient air is subsaturated with respect to ice
  (relative humidity over ice, RHi < 1).  Segments are deactivated once
  their ice water path has dropped to near zero.

Typical usage::

    import numpy as np
    from datetime import datetime, timedelta
    from opendrift.models.contrail import ContrailDrift

    c = ContrailDrift(loglevel=20)
    c.add_readers_from_list(['nwp_forecast.nc'])

    # Seed a contrail along a trans-Atlantic flight path at 11 km
    lons = np.linspace(-10, -60, 200)
    lats = np.linspace(51, 45, 200)
    c.seed_elements(lon=lons, lat=lats, z=11000,
                    time=datetime(2024, 6, 1, 12, 0))

    c.run(duration=timedelta(hours=8), time_step=300,
          time_step_output=1800)
    c.plot(fast=True)
"""

import numpy as np
import logging

logger = logging.getLogger(__name__)

from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements import LagrangianArray
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED

# Thermodynamic constants
_L_SUB = 2.83e6    # Latent heat of sublimation [J kg-1]
_R_V   = 461.5     # Specific gas constant for water vapour [J kg-1 K-1]
_T0    = 273.16    # Triple-point temperature [K]
_E0    = 611.73    # Saturation vapour pressure at triple point [Pa]
_R_D   = 287.05    # Specific gas constant for dry air [J kg-1 K-1]

# SAC mixing-line constants (kerosene / Jet-A fuel)
_CP_AIR  = 1004.0   # Specific heat of air at constant pressure [J kg-1 K-1]
_EPS     = 0.622    # Ratio of molar masses Mw / Ma [-]
_EI_H2O  = 1.25     # Emission index for water vapour [kg_H2O kg_fuel-1]
_Q_FUEL  = 43.2e6   # Specific combustion heat of kerosene [J kg-1]
_ETA     = 0.30     # Typical overall propulsion efficiency [-]


class ContrailElement(LagrangianArray):
    """Lagrangian element representing a contrail ice-crystal segment."""

    variables = LagrangianArray.add_variables([
        ('ice_water_path', {
            'dtype': np.float32,
            'units': 'kg/m2',
            'description': (
                'Vertically integrated ice water content per unit area '
                'of the contrail segment.  Decreases through sublimation.'),
            'level': CONFIG_LEVEL_ESSENTIAL,
            'default': 1e-4,
        }),
        ('optical_depth', {
            'dtype': np.float32,
            'units': '1',
            'description': (
                'Visible optical depth of the contrail segment. '
                'Scales proportionally with ice_water_path.'),
            'level': CONFIG_LEVEL_BASIC,
            'default': 0.1,
        }),
        ('fall_speed', {
            'dtype': np.float32,
            'units': 'm/s',
            'description': (
                'Terminal fall speed of ice crystals (positive value = '
                'downward motion).  Typical range 0.01–0.10 m/s.'),
            'level': CONFIG_LEVEL_ADVANCED,
            'default': 0.03,
        }),
    ])


class ContrailDrift(OpenDriftSimulation):
    """Lagrangian contrail transport and microphysics model.

    Simulates the drift, spreading and sublimation of aircraft contrail
    segments in the upper troposphere.  Each Lagrangian element (particle)
    represents a short contrail segment characterised by its geographic
    position, altitude, ice water path, optical depth and ice-crystal fall
    speed.

    Environmental forcing
    ---------------------
    The model requires the following CF-convention variable names from
    attached readers (fallback defaults shown):

    * ``x_wind`` / ``y_wind`` — horizontal wind components [m s-1].
    * ``air_temperature`` — ambient air temperature [K] (fallback: 220 K).
    * ``specific_humidity`` — water-vapour specific humidity [kg kg-1]
      (fallback: 0, i.e. completely dry → immediate sublimation).
    * ``horizontal_diffusivity`` — horizontal eddy diffusivity [m² s-1]
      (fallback: 50 m² s-1).

    Configuration
    -------------
    ``contrail:sublimation_timescale`` (default 1800 s)
        Characteristic e-folding time for ice sublimation when the ambient
        air is completely dry (RHi = 0).

    ``contrail:sac_eta`` (default 0.30)
        Overall propulsion efficiency η used in the Schmidt-Appleman
        Criterion (SAC).  Higher efficiency → lower T_SAC → contrails
        less likely.

    ``contrail:sac_ei_h2o`` (default 1.25 kg kg-1)
        Emission index for water vapour (kg of H₂O per kg of fuel burned).

    ``contrail:sac_q_fuel`` (default 43.2 × 10⁶ J kg-1)
        Specific combustion heat of the fuel (kerosene / Jet-A).

    Notes
    -----
    * Seed elements with ``z > 0`` (metres above sea level).  Cruise
      altitudes are typically 8 000–12 000 m.
    * Coastline interaction is disabled by default (contrails fly over land).
    * ``seed:ocean_only`` is set to ``False`` so elements can be seeded
      above land.
    * Elements are deactivated with reason ``'no_contrail'`` when the SAC
      is not satisfied, and with reason ``'sublimated'`` when all ice has
      sublimated.
    """

    ElementType = ContrailElement

    required_variables = {
        'x_wind':                {'fallback': 0},
        'y_wind':                {'fallback': 0},
        'air_temperature':       {'fallback': 220},    # K
        'specific_humidity':     {'fallback': 0},      # kg kg-1
        'horizontal_diffusivity':{'fallback': 50},     # m2 s-1
    }

    status_colors = {
        'active':      'cyan',
        'sublimated':  'lightgray',
        'no_contrail': 'orange',
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._add_config({
            'contrail:sublimation_timescale': {
                'type': 'float',
                'default': 1800.,
                'min': 60.,
                'max': 86400.,
                'units': 's',
                'description': (
                    'Characteristic sublimation time scale when ambient '
                    'relative humidity over ice is zero.  The ice loss rate '
                    'scales linearly with the RHi deficit.'),
                'level': CONFIG_LEVEL_ADVANCED,
            },
            'contrail:sac_eta': {
                'type': 'float',
                'default': _ETA,
                'min': 0.05,
                'max': 0.70,
                'units': '-',
                'description': (
                    'Overall propulsion efficiency η for the Schmidt-Appleman '
                    'Criterion (SAC).  Higher η reduces the mixing-line slope '
                    'G and lowers T_SAC, making contrail formation less likely.'),
                'level': CONFIG_LEVEL_ADVANCED,
            },
            'contrail:sac_ei_h2o': {
                'type': 'float',
                'default': _EI_H2O,
                'min': 0.5,
                'max': 2.0,
                'units': 'kg kg-1',
                'description': (
                    'Water-vapour emission index for the Schmidt-Appleman '
                    'Criterion: kg of H₂O per kg of fuel burned.'),
                'level': CONFIG_LEVEL_ADVANCED,
            },
            'contrail:sac_q_fuel': {
                'type': 'float',
                'default': _Q_FUEL,
                'min': 30.e6,
                'max': 50.e6,
                'units': 'J kg-1',
                'description': (
                    'Specific combustion heat of the fuel used in the '
                    'Schmidt-Appleman Criterion.'),
                'level': CONFIG_LEVEL_ADVANCED,
            },
        })

        # Contrails are atmospheric; no coastline or ocean-only constraints.
        self._set_config_default('general:coastline_action', 'none')
        self._set_config_default('seed:ocean_only', False)
        self._set_config_default('drift:max_speed', 150)   # allow jet-stream

    # ------------------------------------------------------------------
    # Thermodynamics helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _saturation_pressure_over_ice(T):
        """Saturation vapour pressure over ice [Pa] via Clausius–Clapeyron.

        Uses the standard meteorological approximation::

            e_sat = e0 * exp(L_sub / R_v * (1/T0 - 1/T))

        Parameters
        ----------
        T : array-like
            Air temperature [K].

        Returns
        -------
        numpy.ndarray
            Saturation vapour pressure [Pa].
        """
        T = np.asarray(T, dtype=float)
        return _E0 * np.exp(_L_SUB / _R_V * (1.0 / _T0 - 1.0 / np.maximum(T, 150.)))

    @staticmethod
    def _saturation_pressure_over_liquid(T):
        """Saturation vapour pressure over liquid water [Pa] (Magnus formula).

        Parameters
        ----------
        T : array-like
            Air temperature [K].

        Returns
        -------
        numpy.ndarray
            Saturation vapour pressure over liquid [Pa].
        """
        T = np.asarray(T, dtype=float)
        T_C = T - 273.15
        return 611.2 * np.exp(17.67 * T_C / (T_C + 243.5))

    @staticmethod
    def _sac_critical_temperature(G):
        """Critical temperature T_SAC [K] for the Schmidt-Appleman Criterion.

        Finds T_SAC by solving ``d(e_sat_liq)/dT = G`` via Newton–Raphson
        iteration on the Magnus formula.  The solution is the temperature at
        which the exhaust-plume mixing line in the (T, e) diagram is tangent
        to the liquid-water saturation curve.

        Parameters
        ----------
        G : array-like
            Mixing-line slope [Pa K-1].

        Returns
        -------
        numpy.ndarray
            Critical temperature T_SAC [K].
        """
        a, b, c = 17.67, 243.5, 611.2   # Magnus coefficients; e [Pa], T [°C]
        G = np.asarray(G, dtype=float)
        T = np.full_like(G, 240.0)       # initial guess [K]
        for _ in range(12):
            T_C     = T - 273.15
            e       = c * np.exp(a * T_C / (T_C + b))
            de_dT   = e * a * b / (T_C + b) ** 2
            d2e_dT2 = de_dT * (a * b / (T_C + b) ** 2 - 2.0 / (T_C + b))
            step    = np.where(np.abs(d2e_dT2) > 1e-15,
                               (de_dT - G) / d2e_dT2, 0.0)
            T      -= step
            T       = np.clip(T, 195.0, 273.15)
        return T

    def _relative_humidity_over_ice(self):
        """Compute relative humidity over ice (RHi) for all active elements.

        Ambient pressure is estimated from the International Standard
        Atmosphere (ISA) using element altitude ``z``.

        Returns
        -------
        numpy.ndarray
            RHi [-], one value per active element.  Values > 1 indicate
            supersaturation; values < 1 indicate subsaturation.
        """
        T = self.environment.air_temperature      # K
        q = self.environment.specific_humidity    # kg kg-1

        # ISA pressure from altitude [Pa]: valid below ~51 km
        z = np.maximum(self.elements.z, 0.)
        p = 101325. * np.maximum(1. - 2.2557e-5 * z, 0.) ** 5.2559

        # Partial pressure of water vapour
        e = q * p / (0.622 + 0.378 * q)

        e_sat = self._saturation_pressure_over_ice(T)
        return np.where(e_sat > 0., e / e_sat, 0.)

    # ------------------------------------------------------------------
    # Update method and sub-steps
    # ------------------------------------------------------------------

    def update(self):
        """Advance contrail elements by one time step.

        Performs the following sub-steps in order:

        1. Schmidt-Appleman Criterion check: deactivate elements where
           ambient temperature exceeds T_SAC (no contrail can form).
        2. Horizontal advection with upper-level winds.
        3. Horizontal diffusion (random walk).
        4. Gravitational sedimentation (lower altitude by fall speed × dt).
        5. Sublimation: reduce ice water path and optical depth when
           RHi < 1.
        6. Deactivate elements whose ice has fully sublimated.
        """
        self._check_sac()
        self._advect_with_wind()
        self.horizontal_diffusion()
        self._sediment()
        self._sublimate()
        self._deactivate_sublimated()

    def _check_sac(self):
        """Deactivate elements where the Schmidt-Appleman Criterion is not met.

        The SAC states that contrails form only when ambient temperature is
        below the critical temperature T_SAC, i.e. when the exhaust-plume
        mixing line in the (T, e) diagram can reach liquid-water saturation.

        The mixing-line slope G [Pa K-1] depends on aircraft parameters and
        local pressure::

            G = cp_air * p * EI_H2O / (eps * Q_fuel * (1 - eta))

        T_SAC is found by solving de_sat_liq/dT = G via Newton–Raphson
        iteration.

        Elements with T > T_SAC are deactivated with reason ``'no_contrail'``.
        """
        eta    = self.get_config('contrail:sac_eta')
        EI_H2O = self.get_config('contrail:sac_ei_h2o')
        Q      = self.get_config('contrail:sac_q_fuel')

        # ISA pressure at element altitude [Pa]
        z = np.maximum(self.elements.z, 0.)
        p = 101325. * np.maximum(1. - 2.2557e-5 * z, 0.) ** 5.2559

        # Mixing-line slope G [Pa K-1]
        G = _CP_AIR * p * EI_H2O / (_EPS * Q * (1.0 - eta))

        T_SAC = self._sac_critical_temperature(G)
        T     = self.environment.air_temperature

        no_contrail = T > T_SAC
        if np.any(no_contrail):
            logger.debug(
                '%d elements deactivated: SAC not satisfied '
                '(T_mean=%.1f K, T_SAC_mean=%.1f K)' % (
                    int(np.sum(no_contrail)),
                    float(T[no_contrail].mean()),
                    float(T_SAC[no_contrail].mean())))
            self.deactivate_elements(no_contrail, reason='no_contrail')

    def _advect_with_wind(self):
        """Advect contrail segments horizontally with upper-tropospheric wind."""
        self.update_positions(self.environment.x_wind,
                              self.environment.y_wind)

    def _sediment(self):
        """Lower elements by their ice-crystal terminal fall speed."""
        dt = self.time_step.total_seconds()
        self.elements.z -= self.elements.fall_speed * dt
        logger.debug('Sedimentation: elements lowered by %.2f m' %
                     (self.elements.fall_speed.mean() * dt))

    def _sublimate(self):
        """Reduce ice water path in subsaturated ambient air.

        The ice water path decays exponentially with a time scale
        ``tau = sublimation_timescale / max(1 - RHi, epsilon)``.
        At full saturation (RHi ≥ 1) the ice is stable and no sublimation
        occurs.
        """
        RHi    = self._relative_humidity_over_ice()
        deficit = np.maximum(1. - RHi, 0.)   # subsaturation fraction [0, 1]

        tau_sub = self.get_config('contrail:sublimation_timescale')
        dt      = self.time_step.total_seconds()

        # Fractional loss over this time step (linear approximation of exp decay)
        loss = np.minimum(deficit * dt / tau_sub, 1.)

        self.elements.ice_water_path *= (1. - loss)
        self.elements.optical_depth  *= (1. - loss)

        n_sub = np.sum(deficit > 0)
        if n_sub > 0:
            logger.debug('%d elements sublimating (mean RHi deficit %.3f)' %
                         (n_sub, deficit[deficit > 0].mean()))

    def _deactivate_sublimated(self):
        """Deactivate elements whose ice water path has dropped to zero."""
        fully_sublimated = self.elements.ice_water_path < 1e-8
        if np.any(fully_sublimated):
            logger.debug('%d contrail elements fully sublimated' %
                         np.sum(fully_sublimated))
            self.deactivate_elements(fully_sublimated, reason='sublimated')
