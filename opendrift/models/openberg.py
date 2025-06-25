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
# Copyright 2015, 2023, Knut-Frode Dagestad, MET Norway
# Copyright 2024, Lenny Hucher, NERSC, Norway
# Copyright 2023, 2024, 2025 Achref Othmani, NERSC, Norway

"""
This code is initiated from the following reference with posterior modifications. 

Reference:
Keghouche, I., F. Counillon, and L. Bertino (2010), Modeling dynamics and thermodynamics
of icebergs in the Barents Sea from 1987 to 2005, J. Geophys. Res., 115, C12062, doi:10.1029/2010JC006165. 
"""

import logging; logger = logging.getLogger(__name__)
from opendrift.elements import LagrangianArray
from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.config import CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ESSENTIAL
from opendrift.models.physics_methods import PhysicsMethods
from scipy.integrate import solve_ivp
import numpy as np


# Constants
rho_water = 1027       # Density of water (kg/m^3)
rho_air = 1.293        # Density of air (kg/m^3)
rho_ice = 917          # Density of ice (kg/m^3)
rho_iceb = 900         # Density of iceberg (kg/m^3)
g = 9.81               # Acceleration due to gravity in m/s²
omega = 7.2921e-5      # Angular frequency (rad/s)
csi = 1                # Sea ice coefficient of resistance
wave_drag_coef = 0.3   # Wave drag coefficient



class IcebergObj(LagrangianArray):
    """ Extending LagrangianArray with relevant properties for an Iceberg """
    variables = LagrangianArray.add_variables([
        ('sail', {'dtype': np.float32,
                  'units': 'm',
                  'default': 10,
                  'description': 'Sail of iceberg (part above waterline)',
                  'level': CONFIG_LEVEL_ESSENTIAL}),
        ('draft', {'dtype': np.float32,
                   'units': 'm',
                   'default': 90,
                   'description': 'Draft of iceberg (part below waterline)',
                   'level': CONFIG_LEVEL_ESSENTIAL}),
        ('length', {'dtype': np.float32,
                    'units': 'm',
                    'default': 100,
                    'description': 'Length of iceberg',
                    'level': CONFIG_LEVEL_ESSENTIAL}),
        ('width', {'dtype': np.float32,
                   'units': 'm',
                   'default': 30,
                   'description': 'Width of iceberg)',
                   'level': CONFIG_LEVEL_ESSENTIAL}),
        ('weight_coeff', {'dtype': np.float32,     # Relative to the shape of iceberg (e.g. 1 for tabular; 0.3 for pinnacle: It affects the mass only)
                              'units': '1',
                              'default': 1}),
        ('water_drag_coeff', {'dtype': np.float32, # Ocean drag coeff.
                              'units': '1',
                              'default': 0.25}),
        ('wind_drag_coeff', {'dtype': np.float32,  # Wind/Air drag coeff.
                             'units': '1',
                             'default': 0.7}),
        ("iceb_x_velocity", {"dtype": np.float32,  # Iceberg velocity in the x-direction
                             "units": "m/s",
                             "default": 0.0}),
        ("iceb_y_velocity", {"dtype": np.float32,  # Iceberg velocity in the y-direction
                             "units": "m/s",
                             "default": 0.0}),
        ])


# Define the functions needed
def ocean_force(iceb_vel, water_vel, Ao, rho_water, water_drag_coef):
    """ Ocean force
    Args:
        iceb_vel  : Iceberg's velocity at time t
        water_vel : Ocean current velocity
        Ao : Iceberg's area in contact with ocean (length x draft)
        rho_water : Water density
        water_drag_coef : Co is the drag coefficient applied on the iceberg's draft
    """
    vxo, vyo = water_vel[0], water_vel[1]
    x_vel, y_vel = iceb_vel[0], iceb_vel[1]
    rel_water_x_vel = vxo - x_vel
    rel_water_y_vel = vyo - y_vel
    rel_water_norm = np.sqrt(rel_water_x_vel**2 + rel_water_y_vel**2)
    F_ocean_x = (0.5 * rho_water * water_drag_coef * Ao * rel_water_norm * rel_water_x_vel)
    F_ocean_y = (0.5 * rho_water * water_drag_coef * Ao * rel_water_norm * rel_water_y_vel)
    return np.array([F_ocean_x, F_ocean_y])


def wind_force(iceb_vel, wind_vel, Aa, wind_drag_coef):
    """ Wind force
    Args:
        iceb_vel : Iceberg's velocity at time t
        wind_vel : Wind velocity
        Aa : Iceberg's area in contact with wind (length x sail)
        wind_drag_coef : Ca is the drag coefficient applied on the iceberg's sail
    """
    vxa, vya = wind_vel[0], wind_vel[1]
    x_vel, y_vel = iceb_vel[0], iceb_vel[1]
    rel_wind_x_vel = vxa - x_vel
    rel_wind_y_vel = vya - y_vel
    rel_wind_norm = np.sqrt(rel_wind_x_vel**2 + rel_wind_y_vel**2)
    F_wind_x = 0.5 * rho_air * wind_drag_coef * Aa * rel_wind_norm * rel_wind_x_vel
    F_wind_y = 0.5 * rho_air * wind_drag_coef * Aa * rel_wind_norm * rel_wind_y_vel
    return np.array([F_wind_x, F_wind_y])


def wave_radiation_force(rho_water, wave_height, wave_direction, iceb_length):
    """ Wave radiation force
    Args:
        rho_water : Water density
        wave_height    : Wave significant height
        wave_direction : Wave direction
        iceb_length    : Iceberg's length
    """
    F_wave_x = (0.5 * rho_water * wave_drag_coef * g * iceb_length * (wave_height / 2) ** 2 * np.sin(np.deg2rad(wave_direction)))
    F_wave_y = (0.5 * rho_water * wave_drag_coef * g * iceb_length * (wave_height / 2) ** 2 * np.cos(np.deg2rad(wave_direction)))
    return np.array([F_wave_x, F_wave_y])


def advect_iceberg_no_acc(f, water_vel, wind_vel):
    """ Advect iceberg without acceleration
    Args:
        f : Wind drift factor
        water_vel : Ocean current velocity
        wind_vel  : Wind velocity
    Returns:
        Iceberg's velocity without acceleration
    """
    vxo, vyo = water_vel[0], water_vel[1]
    vxa, vya = wind_vel[0], wind_vel[1]
    no_acc_vel_x = (1 - f) * vxo + f * vxa
    no_acc_vel_y = (1 - f) * vyo + f * vya
    V = np.array([no_acc_vel_x, no_acc_vel_y])
    if not np.isfinite(V).all():
        logger.error("Infinite value in iceberg's velocity without acceleration: Please check the wind drift factor f ")
    return V


def sea_ice_force(iceb_vel, sea_ice_conc, Ai, sea_ice_vel, sum_force):
    """ Sea ice force
    Args:
        iceb_vel : Iceberg velocity at time t
        sea_ice_conc : Sea ice concentration
        Ai : Iceberg's area in contact with ice (sea_ice_thickness x length) # (Alternatively: Test half length and half width)
        sea_ice_vel : Sea ice velocity
        sum_force : Effect of all other forces exerted on the iceberg (apart from the sea ice force)
    """
    ice_x, ice_y = sea_ice_vel
    x_vel, y_vel = iceb_vel[0], iceb_vel[1]
    diff_vel = np.sqrt((ice_x - x_vel) ** 2 + (ice_y - y_vel) ** 2)
    force_x, force_y = sum_force
    F_ice_x = np.zeros_like(x_vel)
    F_ice_y = np.zeros_like(y_vel)
    F_ice_x = (0.5 * (rho_ice * csi * Ai) * diff_vel * (ice_x - x_vel))
    F_ice_y = (0.5 * (rho_ice * csi * Ai) * diff_vel * (ice_y - y_vel))
    F_ice_x[sea_ice_conc <= 0.15] = 0
    F_ice_y[sea_ice_conc <= 0.15] = 0
    F_ice_x[sea_ice_conc >= 0.9] = -force_x[sea_ice_conc >= 0.9]
    F_ice_y[sea_ice_conc >= 0.9] = -force_y[sea_ice_conc >= 0.9]
    return np.array([F_ice_x, F_ice_y])


def coriolis_force(iceb_vel, mass, lat):
    """ Coriolis force
    Args:
        iceb_vel : Iceberg velocity at time t
        mass: Mass of the iceberg
        lat : Latitude of the iceberg's location in degrees
    """
    f = 2 * omega * np.sin(np.radians(lat))
    assert len(iceb_vel) == 2
    x_vel, y_vel = iceb_vel[0], iceb_vel[1]
    F_cor_x = mass * f * y_vel
    F_cor_y = -mass * f * x_vel
    return np.array([F_cor_x, F_cor_y])


def sea_surface_slope_force(sea_slope_x, sea_slope_y, mass):
    """ This functions assumes you provide the sea surface slope from an external file """
    F_sea_slope_x = -mass * g * sea_slope_x
    F_sea_slope_y = mass * g * sea_slope_y
    return np.array([F_sea_slope_x, F_sea_slope_y])


def melwav(iceb_length, iceb_width, x_wind, y_wind, sst, sea_ice_conc, dt):
    """ Update the iceberg's dimensions (length and width) due to wave erosion
    Args:
        iceb_length : Iceberg's length
        iceb_width : Iceberg's width
        x_wind : Wind speed in the x-direction
        y_wind : Wind speed in the y-direction
        sst : Sea surface temperature
        sea_ice_conc : Sea ice concentration
        dt : Timestep of the simulation
    """
    Ss = -5 + np.sqrt(32 + 2 * np.sqrt(x_wind**2 + y_wind**2))
    Vsst = (1 / 6.0) * (sst + 2) * Ss
    Vwe = Vsst * 0.5 * (1 + np.cos(np.pi * sea_ice_conc**3)) / 86400
    new_iceb_length = np.zeros_like(iceb_length)
    new_iceb_width = np.zeros_like(iceb_width)
    new_iceb_length[iceb_length != 0] = (iceb_length[iceb_length != 0] - Vwe[iceb_length != 0] * dt)
    new_iceb_width[iceb_length != 0] = (iceb_width[iceb_length != 0] / iceb_length[iceb_length != 0] * new_iceb_length[iceb_length != 0])
    new_iceb_length[new_iceb_length < 0] = 0
    new_iceb_width[new_iceb_width < 0] = 0
    return new_iceb_length, new_iceb_width


def mellat(iceb_length, iceb_width, tempib, salnib, dt): 
    """ Update the iceberg's dimensions (length and width) due to lateral melting
    Args:
        iceb_length : Iceberg's length
        iceb_width  : Iceberg's width
        tempib : Water temperature
        salnib : Water salinity
        dt : Timestep of the simulation
    """
    TfS = -0.036 - 0.0499 * salnib - 0.000112 * salnib**2
    Tfp = TfS * np.exp(-0.19 * (tempib - TfS))
    deltaT = tempib - Tfp
    deltaT = np.concatenate([2.78 * deltaT, 0.47 * deltaT**2], axis=0)
    sumVb = np.nansum(deltaT, axis=0)
    dx = sumVb / 365 / 86400 * dt
    
    new_iceb_length = np.zeros_like(iceb_length)
    new_iceb_width = np.zeros_like(iceb_width)
    new_iceb_length[iceb_length != 0] = (iceb_length[iceb_length != 0] - 2 * dx[iceb_length != 0])
    new_iceb_width[iceb_length != 0] = (iceb_width[iceb_length != 0] / iceb_length[iceb_length != 0] * new_iceb_length[iceb_length != 0])
    new_iceb_length[new_iceb_length < 0] = 0
    new_iceb_width[new_iceb_width < 0] = 0
    return new_iceb_length, new_iceb_width


def melbas(iceb_draft, iceb_sail, iceb_length, salnib, tempib, x_water_vel, y_water_vel, x_iceb_vel, y_iceb_vel, dt):
    """ Update the iceberg's dimensions (draft and sail) due to forced convection """
    # Temperature at the base layer of the iceberg
    absv = np.sqrt(((x_water_vel - x_iceb_vel) ** 2 + (y_water_vel - y_iceb_vel) ** 2))
    TfS = -0.036 - 0.0499 * salnib - 0.000112 * salnib**2
    Tfp = TfS * 2.71828 ** (-0.19 * (tempib - TfS))
    deltat = tempib - Tfp
    Vf = 0.58 * absv**0.8 * deltat / (iceb_length**0.2)
    Vf = Vf / 86400  # conversion to m/s
    # Update the draft
    new_iceb_draft = np.zeros_like(iceb_draft)
    new_iceb_draft[iceb_draft != 0] = (abs(iceb_draft[iceb_draft != 0]) - Vf[iceb_draft != 0] * dt)
    # Melt at base of the iceberg
    new_iceb_draft[iceb_draft < 0] = 0
    return new_iceb_draft, iceb_sail


class OpenBerg(OpenDriftSimulation):

    ElementType = IcebergObj

    required_variables = {
        "x_sea_water_velocity": {"fallback": None, "profiles": True},
        "y_sea_water_velocity": {"fallback": None, "profiles": True},
        "sea_floor_depth_below_sea_level": {"fallback": 10000},
        'sea_surface_height': {'fallback': 0},
        "sea_surface_x_slope": {"fallback": 0},
        "sea_surface_y_slope": {"fallback": 0},
        "x_wind": {"fallback": None, "important": True},
        "y_wind": {"fallback": None, "important": True},
        "sea_surface_wave_significant_height": {"fallback": 0},
        "sea_surface_wave_from_direction": {"fallback": 0},
        "sea_surface_wave_stokes_drift_x_velocity": {"fallback": 0},
        "sea_surface_wave_stokes_drift_y_velocity": {"fallback": 0},
        "sea_water_temperature": {"fallback": 2, "profiles": True},
        "sea_water_salinity": {"fallback": 35, "profiles": True},
        "sea_ice_area_fraction": {"fallback": 0},
        "sea_ice_thickness": {"fallback": 0},
        "sea_ice_x_velocity": {"fallback": 0, "important": False},
        "sea_ice_y_velocity": {"fallback": 0, "important": False},
        "land_binary_mask": {"fallback": None},
    }


    def get_profile_masked(self, variable):
        """
        Apply a mask to extract data from the surface down to the iceberg's draft.
        """
        draft = self.elements.draft
        profile = self.environment_profiles[variable]
        z = self.environment_profiles["z"] 
        if profile.ndim == 1:
            profile = profile[np.newaxis, :]
        if z is None or (len(z) == 1 and z[0] is None):
            z = np.zeros(profile.shape[0])
        z = np.atleast_1d(z)
        if z.ndim == 1:
            z = z[:, np.newaxis] 
        draft = np.atleast_1d(draft)
        mask = draft[np.newaxis, :] < -z 
        if mask.shape[0] > 1:
            mask[np.argmax(mask, axis=0), np.arange(mask.shape[1])] = False
        assert profile.shape == mask.shape, f"Incompatible shapes: profile {profile.shape}, mask {mask.shape}"

        return np.ma.masked_array(profile, mask, fill_value=np.nan)


    def get_basal_env(self, variable):
        """ Get the basal layer of the variable for the icebergs """
        profile = self.get_profile_masked(variable)
        last = np.argmin(np.logical_not(profile.mask), axis=0) - 1
        return profile[last, np.arange(profile.shape[1])]


    # Configuration
    def __init__(self, *args, **kwargs):
        super(OpenBerg, self).__init__(*args, **kwargs)
        self._add_config({
            'drift:wave_rad':{
                'type': 'bool',
                'default': True,
                'description': 'If True, wave radiation force is added',
                'level': CONFIG_LEVEL_BASIC
            },
            'drift:stokes_drift':{
                'type': 'bool',
                'default': False,
                'description': 'If True, stokes drift force is added',
                'level': CONFIG_LEVEL_BASIC
            },
            'drift:coriolis':{
                'type': 'bool',
                'default': True,
                'description': 'If True, coriolis force is added',
                'level': CONFIG_LEVEL_BASIC,
            },
            'drift:sea_surface_slope':{
            'type': 'bool',
            'default': False,
            'description': 'If True, sea surface slope force is added',
            'level': CONFIG_LEVEL_BASIC,
            },
            'drift:vertical_profile':{
                'type': 'bool',
                'default': False,
                'description': 'If True, depth integrated currents are applied',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:grounding':{
                'type': 'bool',
                'default': True,
                'description': 'If True, grounding is enabled',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:roll_over':{
                'type': 'bool',
                'default': True,
                'description': 'If True, roll over is enabled',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:melting':{
                'type': 'bool',
                'default': False,
                'description': 'If True, melting is enabled',
                'level': CONFIG_LEVEL_BASIC
            },
            'melting:wave':{
                'type': 'bool',
                'default': True,
                'description': 'If True, melting due to wave erosion is enabled',
                'level': CONFIG_LEVEL_BASIC
            },            
            'melting:lateral':{
                'type': 'bool',
                'default': True,
                'description': 'If True, lateral melting is enabled',
                'level': CONFIG_LEVEL_BASIC
            },  
            'melting:basal':{
                'type': 'bool',
                'default': True,
                'description': 'If True, basal melting is enabled',
                'level': CONFIG_LEVEL_BASIC
            },  
        })


    def advect_iceberg(self):
        sail = self.elements.sail
        draft = self.elements.draft
        length = self.elements.length
        width = self.elements.width
        weight_coeff = self.elements.weight_coeff
        lat = self.elements.lat
        water_drag_coeff = self.elements.water_drag_coeff
        wind_drag_coeff = self.elements.wind_drag_coeff

        T = self.environment.sea_water_temperature
        S = self.environment.sea_water_salinity
        rho_water = PhysicsMethods.sea_water_density(T, S)
        sea_slope_x = self.environment.sea_surface_x_slope
        sea_slope_y = self.environment.sea_surface_y_slope
        sea_surface_height= self.environment.sea_surface_height
        wave_height = self.environment.sea_surface_wave_significant_height
        wave_direction = self.environment.sea_surface_wave_from_direction
        sea_ice_thickness = self.environment.sea_ice_thickness
        sea_ice_conc = self.environment.sea_ice_area_fraction
        water_depth = self.environment.sea_floor_depth_below_sea_level

        Ao = abs(draft) * length # (Alternatively: Ao = weight_coeff * length * width)
        Aa = sail * length
        Ai = sea_ice_thickness * length
        mass = width * (Aa + Ao) * rho_iceb * weight_coeff
        k = (rho_air * wind_drag_coeff * Aa / (rho_water * water_drag_coeff * Ao))
        f = np.sqrt(k) / (1 + np.sqrt(k)) # (f is the wind drift factor, only used in the no acceleration model)

        wave_rad = self.get_config('drift:wave_rad')
        stokes_drift = self.get_config('drift:stokes_drift')
        coriolis = self.get_config('drift:coriolis')
        grounding = self.get_config('processes:grounding')
        sea_surface_slope = self.get_config('drift:sea_surface_slope')
        
                
        if self.get_config('drift:vertical_profile') is False:
            logger.debug("Advection with surface currents")
            water_vel = np.array([self.environment.x_sea_water_velocity + (int(stokes_drift) * self.environment.sea_surface_wave_stokes_drift_x_velocity),
                                  self.environment.y_sea_water_velocity + (int(stokes_drift) * self.environment.sea_surface_wave_stokes_drift_y_velocity)])
        else:
            logger.debug("Advection with depth integrated currents")
            uprof = self.get_profile_masked("x_sea_water_velocity")
            vprof = self.get_profile_masked("y_sea_water_velocity")
            z = self.environment_profiles["z"]
            thickness = -(z[1:] - z[:-1]).reshape((-1, 1)).astype(float)
            mask = uprof.mask
            uprof_mean_inter = (uprof[1:] + uprof[:-1]) / 2
            vprof_mean_inter = (vprof[1:] + vprof[:-1]) / 2
            mask = mask[:-1]
            thickness_reshaped = np.tile(thickness, (1, mask.shape[1]))
            thickness_reshaped = np.where(mask, np.nan, thickness_reshaped)
            umean = np.nansum(thickness_reshaped * uprof_mean_inter, axis=0) / np.nansum(thickness_reshaped, axis=0)
            vmean = np.nansum(thickness_reshaped * vprof_mean_inter, axis=0) / np.nansum(thickness_reshaped, axis=0)
            water_vel = np.array([umean, vmean])

        wind_vel = np.array([self.environment.x_wind, self.environment.y_wind])
        sea_ice_vel = np.array([self.environment.sea_ice_x_velocity, self.environment.sea_ice_y_velocity])


        def dynamic(t,iceb_vel, water_vel, wind_vel, wave_height, wave_direction, Ao,
                    Aa, rho_water, water_drag_coef, wind_drag_coef, iceb_length, mass,lat, sea_slope_x, sea_slope_y):
            """ Function required by solve_ivp. The t and iceb_vel parameters are required by solve_ivp, shouldn't be deleted """
            iceb_vel = iceb_vel.reshape((2, -1))
            # Individual forces
            ocean_force_val = ocean_force(iceb_vel, water_vel, Ao, rho_water, water_drag_coef)
            wind_force_val = wind_force(iceb_vel, wind_vel, Aa, wind_drag_coef)
            wave_radiation_force_val = int(wave_rad) * wave_radiation_force(rho_water, wave_height, wave_direction, iceb_length)
            coriolis_force_val = int(coriolis) * coriolis_force(iceb_vel, mass, lat)
            sea_surface_slope_val = int(sea_surface_slope) * sea_surface_slope_force(sea_slope_x, sea_slope_y, mass)
            
            # Sum of the individual forces
            sum_force = (ocean_force_val+ wind_force_val+ wave_radiation_force_val+ coriolis_force_val+ sea_surface_slope_val)
            
            # Add sea ice force
            sea_ice_force_val = sea_ice_force(iceb_vel, sea_ice_conc, Ai, sea_ice_vel, sum_force)
            sum_force += sea_ice_force_val

            return (sum_force / mass)

        # Running the simulation
        V0 = advect_iceberg_no_acc(f, water_vel, wind_vel)  # Approximation of the solution of the dynamic equation for the iceberg velocity
        V0[:, sea_ice_conc >= 0.9] = sea_ice_vel[:, sea_ice_conc >= 0.9]  # With this criterium, the iceberg moves with the sea ice
        V0 = V0.flatten() # V0 needs to be 1D

        effective_water_depth = water_depth + sea_surface_height
        hwall = draft - effective_water_depth
        grounded = np.logical_and(hwall >= 0, grounding)
        if grounding:
            # Determine which icebergs are grounded
            if np.any(grounded):
                logger.debug(f"Grounding condition: Icebergs grounded = {np.sum(grounded)}, "
                            f"hwall = {np.round(hwall[grounded], 3)} meters")
                # Grounded icebergs stop moving
                self.elements.moving[grounded] = 0
            else:
                logger.debug("No grounded icebergs detected in this timestep")

            # Check for Degrounding regardless of whether grounding occurred now
            degrounding = np.logical_and(self.elements.moving == 0, hwall < 0)
            if np.any(degrounding):
                logger.debug(f"Degrounding condition: Icebergs degrounded = {np.sum(degrounding)}, "
                            f"hwall = {np.round(hwall[degrounding], 3)} meters")
                # Degrounded icebergs start moving again
                self.elements.moving[degrounding] = 1
        else:
            logger.debug("Grounding process disabled in configuration")
        
        sol = solve_ivp(dynamic, [0, self.time_step.total_seconds()], V0,
                        args=(water_vel, wind_vel, wave_height, wave_direction, Ao, Aa, rho_water,
                              water_drag_coeff, wind_drag_coeff, length, mass, lat, sea_slope_x, sea_slope_y),
                              vectorized=True,
                              t_eval=np.array([self.time_step.total_seconds()]))
        V = sol.y.reshape((2, -1))
        Vx, Vy = V[0], V[1]
        Vx[grounded] = 0
        Vy[grounded] = 0
        self.update_positions(Vx, Vy)
        self.elements.iceb_x_velocity, self.elements.iceb_y_velocity = Vx, Vy


    def melt(self):
        """ Enable melting """
        if self.get_config('processes:melting') is False:
            logger.debug('Melting is disabled')
            return
        x_wind = self.environment.x_wind
        y_wind = self.environment.y_wind
        uoib = self.get_basal_env("x_sea_water_velocity")
        voib = self.get_basal_env("y_sea_water_velocity")
        T_profile = self.environment_profiles["sea_water_temperature"]
        S_profile = self.environment_profiles["sea_water_salinity"]
        Tn = self.get_basal_env("sea_water_temperature")
        Sn = self.get_basal_env("sea_water_salinity")
        sea_ice_conc = self.environment.sea_ice_area_fraction

        # Wave melting
        if self.get_config('melting:wave'):
            self.elements.length, self.elements.width = melwav(self.elements.length, self.elements.width, x_wind, y_wind, T_profile[0], sea_ice_conc, self.time_step.total_seconds())
        # Lateral melting
        if self.get_config('melting:lateral'):
            self.elements.length, self.elements.width = mellat(self.elements.length, self.elements.width, T_profile, S_profile, self.time_step.total_seconds())
        # Basal melting
        if self.get_config('melting:basal'):
            self.elements.draft, self.elements.sail = melbas(self.elements.draft, self.elements.sail, self.elements.length, Sn, Tn, uoib, voib, self.elements.iceb_x_velocity, self.elements.iceb_y_velocity, self.time_step.total_seconds())
        
        # Deactivate icebergs less than 1 meter
        self.deactivate_elements(self.elements.draft < 1, "Iceberg melted")
        self.deactivate_elements(self.elements.length < 1, "Iceberg melted")
        self.deactivate_elements(self.elements.width < 1, "Iceberg melted")
        self.deactivate_elements(self.elements.sail < 1, "Iceberg melted")


    def roll_over(self):
        """ Iceberg's stability criterium """
        if self.get_config('processes:roll_over') is False:
            logger.debug('Rollover is disabled')
            return
        T = self.environment.sea_water_temperature
        S = self.environment.sea_water_salinity
        rho_water = PhysicsMethods.sea_water_density(T, S)
        L = self.elements.length
        W = self.elements.width
        H = self.elements.draft + self.elements.sail
        alpha = rho_iceb / rho_water
        crit = np.sqrt(6 * alpha * (1 - alpha))
        W, L = np.min([L, W], axis=0), np.max([L, W], axis=0)
        mask = (W / H) < crit
        if any(mask):
            logger.debug(f"Rolling over : {np.sum(mask)} icebergs ...")
            nL, nW, nH = (np.max([L[mask], H[mask]], axis=0), np.min([L[mask], H[mask]], axis=0), W[mask])
            L[mask], W[mask], H[mask] = nL, nW, nH
        depthib = H * alpha
        sailib = H - depthib
        self.elements.length = L
        self.elements.width = W
        self.elements.sail = sailib
        self.elements.draft = depthib

    def update(self):
        """ Update positions and properties of particles """
        self.roll_over()
        self.melt()
        self.advect_iceberg()
