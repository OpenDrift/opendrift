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
# Copyright 2023, Achref Othmani & 2024, Lenny Hucher, NERSC, Norway

import numpy as np
import logging; logger = logging.getLogger(__name__)
from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from opendrift.models.physics_methods import PhysicsMethods
from scipy.integrate import solve_ivp



# Constants
rho_water = 1027
rho_air = 1.293
rho_ice = 917
rho_iceb = 900
g = 9.81
csi = 1  # Sea ice coefficient of resistance
wave_drag_coef = 0.3


class IcebergObj(Lagrangian3DArray):
    """ Extending Lagrangian3DArray with relevant properties for an Iceberg """

    variables = Lagrangian3DArray.add_variables([
        ('sail', {'dtype': np.float32,	           # Sail of iceberg (part above waterline )
                               'units': 'm',
                               'default': 10}),
        ('draft', {'dtype': np.float32,	           # Draft of iceberg (part below waterline)
                               'units': 'm',
                               'default': 90}),
        ('length', {'dtype': np.float32,	       # length of iceberg 
                               'units': 'm',
                               'default': 100}),
        ('width', {'dtype': np.float32,		       # width of iceberg 
                               'units': 'm',
                               'default': 30}),
        ('weight_coeff', {'dtype': np.float32,     # Relative to the shape of iceberg (e.g. 1 for tabular; 0.3 for pinnacle)
                              'units': '1',
                              'default': 1}),
        ('water_drag_coeff', {'dtype': np.float32, # Ocean drag coeff.
                              'units': '1',
                              'default': 0.25}),
        ('wind_drag_coeff', {'dtype': np.float32,  # Wind drag coeff.
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
def water_drag(iceb_vel, water_vel, Ao, rho_water, water_drag_coef):
    """ Ocean force
    Args:
        iceb_vel  : Iceberg velocity at time t
        water_vel : Ocean current velocity
        Ao : Iceberg's area exposed to the ocean current(length x draft)
        rho_water : Water density
        water_drag_coef : Co is the drag coefficient applied on the iceberg's draft
    Returns:
        Ocean force
    """
    vxo, vyo = water_vel[0], water_vel[1]    # water velocity
    x_vel, y_vel = iceb_vel[0], iceb_vel[1]  # iceberg velocity
    rel_water_x_vel = vxo - x_vel
    rel_water_y_vel = vyo - y_vel
    rel_water_norm = np.sqrt(rel_water_x_vel**2 + rel_water_y_vel**2)
    F_ocean_x = (0.5 * rho_water * water_drag_coef * Ao * rel_water_norm * rel_water_x_vel)
    F_ocean_y = (0.5 * rho_water * water_drag_coef * Ao * rel_water_norm * rel_water_y_vel)
    return np.array([F_ocean_x, F_ocean_y])


def wind_drag(iceb_vel, wind_vel, Aa, wind_drag_coef):
    """ Wind force
    Args:
        iceb_vel : Iceberg velocity at time t
        wind_vel : Wind velocity
        Aa : Iceberg's area exposed to the wind (length x sail)
        rho_air : Air density
        wind_drag_coef : Ca is the drag coefficient applied on the iceberg's sail
    Returns:
        Wind force
    """
    vxa, vya = wind_vel[0], wind_vel[1] # wind velocity
    x_vel, y_vel = iceb_vel[0], iceb_vel[1] # iceberg velocity
    rel_wind_x_vel = vxa - x_vel
    rel_wind_y_vel = vya - y_vel
    rel_wind_norm = np.sqrt(rel_wind_x_vel**2 + rel_wind_y_vel**2)
    F_wind_x = 0.5 * rho_air * wind_drag_coef * Aa * rel_wind_norm * rel_wind_x_vel
    F_wind_y = 0.5 * rho_air * wind_drag_coef * Aa * rel_wind_norm * rel_wind_y_vel
    return np.array([F_wind_x, F_wind_y])


def wave_radiation_force(rho_water, wave_height, wave_direction, iceb_length):
    """ Wave radiation force
    Args:
        iceb_vel  : Iceberg velocity at time t
        rho_water : Water density
        wave_height    : Wave significant height
        wave_direction : Wave direction
        iceb_length    : Iceberg's length
    Returns:
        Wave radiation force
    """
    F_wave_x = (0.5 * rho_water * wave_drag_coef * g * iceb_length * (wave_height / 2) ** 2 * np.sin(np.deg2rad((wave_direction + 180) % 360)))
    F_wave_y = (0.5 * rho_water * wave_drag_coef * g * iceb_length * (wave_height / 2) ** 2 * np.cos(np.deg2rad((wave_direction + 180) % 360)))
    return np.array([F_wave_x, F_wave_y])


def advect_iceberg_no_acc(f, water_vel, wind_vel):
    """ Advect iceberg without acceleration
    Args:
        f : Wind drift factor
        water_vel : Ocean current velocity
        wind_vel  : Wind velocity

    Returns:
        Iceberg velocity without acceleration
    """
    vxo, vyo = water_vel[0], water_vel[1]
    vxa, vya = wind_vel[0], wind_vel[1]
    no_acc_model_x = (1 - f) * vxo + f * vxa
    no_acc_model_y = (1 - f) * vyo + f * vya
    V = np.array([no_acc_model_x, no_acc_model_y])
    if not np.isfinite(V).all():
        logger.error("Infinite value in iceberg's velocity without acceleration: Please check the wind drift factor f ")
    return V


def sea_ice_force(iceb_vel, sea_ice_conc, sea_ice_thickness, sea_ice_vel, iceb_width, sum_force):
    """ Sea ice force
    Args:
        iceb_vel : Iceberg velocity at time t
        sea_ice_conc : Sea ice concentration [%]
        sea_ice_thickness : Sea ice thickness [m]
        sea_ice_vel : Sea ice velocity [m/s]
        rho_ice : Sea ice density
        iceb_width : Iceberg width
        sum_force : Effect of all other forces exerted on the iceberg (apart from the sea ice force)
    Returns:
        Sea ice force
    """
    ice_x, ice_y = sea_ice_vel
    x_vel, y_vel = iceb_vel[0], iceb_vel[1]
    diff_vel = np.sqrt((ice_x - x_vel) ** 2 + (ice_y - y_vel) ** 2)
    force_x, force_y = sum_force
    F_ice_x = np.zeros_like(x_vel)
    F_ice_y = np.zeros_like(y_vel)
    F_ice_x = (0.5 * (rho_ice * csi * sea_ice_thickness * iceb_width) * diff_vel * (ice_x - x_vel))
    F_ice_y = (0.5 * (rho_ice * csi * sea_ice_thickness * iceb_width) * diff_vel * (ice_y - y_vel))
    F_ice_x[sea_ice_conc <= 0.15] = 0
    F_ice_y[sea_ice_conc <= 0.15] = 0
    F_ice_x[sea_ice_conc >= 0.9] = -force_x[sea_ice_conc >= 0.9]
    F_ice_y[sea_ice_conc >= 0.9] = -force_y[sea_ice_conc >= 0.9]
    return np.array([F_ice_x, F_ice_y])


def melwav(iceb_length, iceb_width, x_wind, y_wind, sst, sea_ice_conc, dt):
    """ Update the iceberg's dimensions (length and width) due to wave erosion
    Parameterization from Keghouche et al. 2009    
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
    Parameterization from Keghouche et al. 2009 
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
    dx = sumVb / 365 / 86400 * dt  # Convert sumVb from (meter/year) to (meter per second)
    
    new_iceb_length = np.zeros_like(iceb_length)
    new_iceb_width = np.zeros_like(iceb_width)
    new_iceb_length[iceb_length != 0] = (iceb_length[iceb_length != 0] - 2 * dx[iceb_length != 0])
    new_iceb_width[iceb_length != 0] = (iceb_width[iceb_length != 0] / iceb_length[iceb_length != 0] * new_iceb_length[iceb_length != 0])
    new_iceb_length[new_iceb_length < 0] = 0
    new_iceb_width[new_iceb_width < 0] = 0
    return new_iceb_length, new_iceb_width

def melbas(iceb_draft, iceb_sail, iceb_length, salnib, tempib, x_water_vel, y_water_vel, x_iceb_vel, y_iceb_vel, dt):
    """ Update the iceberg's dimensions (draft and sail) due to forced convection
    Parameterization from Keghouche et al. 2009 
    """
    # Temperature at the base layer of the iceberg
    absv = np.sqrt(((x_water_vel - x_iceb_vel) ** 2 + (y_water_vel - y_iceb_vel) ** 2))
    TfS = -0.036 - 0.0499 * salnib - 0.000112 * salnib**2
    Tfp = TfS * 2.71828 ** (-0.19 * (tempib - TfS))
    deltat = tempib - Tfp
    Vf = 0.58 * absv**0.8 * deltat / (iceb_length**0.2)
    Vf = Vf / 86400  # converted to m/s
    # Update the draft
    new_iceb_draft = np.zeros_like(iceb_draft)
    new_iceb_draft[iceb_draft != 0] = (abs(iceb_draft[iceb_draft != 0]) - Vf[iceb_draft != 0] * dt)
    # Melt at base of the iceberg
    new_iceb_draft[iceb_draft < 0] = 0
    return new_iceb_draft, iceb_sail


class OpenBerg(OceanDrift):

    ElementType = IcebergObj

    required_variables = {
        "x_sea_water_velocity": {"fallback": None, "profiles": True},
        "y_sea_water_velocity": {"fallback": None, "profiles": True},
        "sea_floor_depth_below_sea_level": {"fallback": 10000},
        "x_wind": {"fallback": 0, "important": False},
        "y_wind": {"fallback": 0, "important": False},
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
        """ Apply a mask to extract data from the surface down to the iceberg's draft """
        draft = self.elements.draft
        profile = self.environment_profiles[variable]
        z = self.environment_profiles["z"]
        if len(z) == 1:
            if z == np.array([None]):  # Convert None values to zeros for surface current
                z = np.zeros_like(z)
        mask = draft[:, np.newaxis] < -z
        mask = mask.T
        mask[np.argmax(mask, axis=0), np.arange(mask.shape[1])] = False
        return np.ma.masked_array(profile, mask, fill_value=np.nan)
    
    def get_basal_env(self, variable):
        """ Get the basal layer of the variable for the icebergs """
        profile = self.get_profile_masked(variable)
        last = np.argmin(np.logical_not(profile.mask), axis=0) - 1
        return profile[last, np.arange(profile.shape[1])]


    # Configuration
    def __init__(
        self,
        add_stokes_drift: bool = True,
        wave_rad: bool = True,
        grounding: bool = False,
        vertical_profile: bool = False,
        melting: bool = False,
        choose_melting: dict[bool] = {"wave": True, "lateral": True, "basal": True},
        *args,
        **kwargs,
    ):

        # The constructor of parent class must always be called
        # to perform some necessary common initialisation tasks:
        super(OpenBerg, self).__init__(*args, **kwargs)
        self.wave_rad = wave_rad
        self.add_stokes_drift = add_stokes_drift
        self.grounding = grounding
        self.vertical_profile = vertical_profile
        self.melting = (melting)   #Activate the melting
        self.choose_melting = (choose_melting)   #Choose how it melts

    def advect_iceberg(self, rho_water, stokes_drift=True, wave_rad=True, grounding=False, vertical_profile=False):
        draft = self.elements.draft
        length = self.elements.length
        Ao = abs(draft) * length  ### Area_wet
        Aa = self.elements.sail * length  ### Area_dry
        mass = self.elements.width * (Aa + Ao) * rho_iceb
        k = (rho_air * self.elements.wind_drag_coeff * Aa / (rho_water * self.elements.water_drag_coeff * Ao))
        f = np.sqrt(k) / (1 + np.sqrt(k))

        if not vertical_profile:
            water_vel = np.array([self.environment.x_sea_water_velocity + int(stokes_drift) * self.environment.sea_surface_wave_stokes_drift_x_velocity,
                                self.environment.y_sea_water_velocity + int(stokes_drift) * self.environment.sea_surface_wave_stokes_drift_y_velocity])
        else:
            logger.info("Calculating the Depth Integrated Currents ...")
            uprof = self.get_profile_masked("x_sea_water_velocity")
            vprof = self.get_profile_masked("y_sea_water_velocity")

            z = self.environment_profiles["z"]
            thickness = -(z[1:] - z[:-1]).reshape((-1, 1))
            mask = uprof.mask
            uprof_mean_inter = (uprof[1:] + uprof[:-1]) / 2
            vprof_mean_inter = (vprof[1:] + vprof[:-1]) / 2
            mask = mask[:-1]
            thickness_reshaped = np.tile(thickness, (1, mask.shape[1]))
            thickness_reshaped[mask] = np.nan
            umean = np.nansum(thickness_reshaped * uprof_mean_inter, axis=0) / np.nansum(thickness_reshaped, axis=0)
            vmean = np.nansum(thickness_reshaped * vprof_mean_inter, axis=0) / np.nansum(thickness_reshaped, axis=0)
            water_vel = np.array([umean, vmean])

        water_depth = self.environment.sea_floor_depth_below_sea_level
        wind_vel = np.array([self.environment.x_wind, self.environment.y_wind])
        wave_height = self.environment.sea_surface_wave_significant_height
        wave_direction = self.environment.sea_surface_wave_from_direction
        sea_ice_conc = self.environment.sea_ice_area_fraction
        sea_ice_thickness = self.environment.sea_ice_thickness
        sea_ice_vel = np.array([self.environment.sea_ice_x_velocity, self.environment.sea_ice_y_velocity])

        def dynamic(t,iceb_vel, water_vel, wind_vel, wave_height, wave_direction, Ao,
                    Aa, rho_water, water_drag_coef, wind_drag_coef, iceb_length, mass):
            """ Function required by solve_ivp. The t and iceb_vel parameters are required by solve_ivp, shouldn't be deleted """
            sum_force = (water_drag(iceb_vel, water_vel, Ao, rho_water, water_drag_coef)
                         + wind_drag(iceb_vel,wind_vel, Aa, wind_drag_coef)
                         + int(wave_rad) * wave_radiation_force(rho_water, wave_height, wave_direction, iceb_length))
            
            sum_force = sum_force + sea_ice_force(iceb_vel, sea_ice_conc, sea_ice_thickness, sea_ice_vel, self.elements.width, sum_force)
            return 1 / mass * sum_force
        
        V0 = advect_iceberg_no_acc(f, water_vel, wind_vel)  # Approximation of the solution of the dynamic equation for the iceberg velocity
        V0[:, sea_ice_conc >= 0.9] = sea_ice_vel[:, sea_ice_conc >= 0.9]  # With this criterium, the iceberg moves with the sea ice
        V0 = V0.flatten()
        hwall = draft - water_depth
        grounded = np.logical_and(hwall >= 0, grounding)
        if any(grounded) and grounding:
            logger.info(f"Grounding condition : Icebergs grounded = {len(hwall[hwall>0])}, hwall={np.round(hwall[hwall>0],3)} meters")
        
        sol = solve_ivp(dynamic, [0, self.time_step.total_seconds()], V0,
                        args=(water_vel, wind_vel, wave_height, wave_direction,
                              Ao, Aa, rho_water,
                              self.elements.water_drag_coeff,
                              self.elements.wind_drag_coeff,
                              self.elements.length,
                              mass),
                              vectorized=True,
                              t_eval=np.array([self.time_step.total_seconds()]),
                              )

        V = sol.y.reshape((2, -1))
        Vx, Vy = V[0], V[1]
        Vx[grounded] = 0
        Vy[grounded] = 0
        self.update_positions(Vx, Vy)
        self.elements.iceb_x_velocity, self.elements.iceb_y_velocity = Vx, Vy

    def melt(self):
        """ Activate melting """
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
        if self.choose_melting["wave"]:
            self.elements.length, self.elements.width = melwav(self.elements.length, self.elements.width, x_wind, y_wind, T_profile[0], sea_ice_conc, self.time_step.total_seconds())
        # Lateral melting
        if self.choose_melting["lateral"]:
            self.elements.length, self.elements.width = mellat(self.elements.length, self.elements.width, T_profile, S_profile, self.time_step.total_seconds())
        # Basal melting
        if self.choose_melting["basal"]:
            self.elements.draft, self.elements.sail = melbas(self.elements.draft, self.elements.sail, self.elements.length, Sn, Tn, uoib, voib, self.elements.iceb_x_velocity, self.elements.iceb_y_velocity, self.time_step.total_seconds())
        
        # Deactivate elements less than 1 meter
        self.deactivate_elements(self.elements.draft < 1, "Iceberg melted")
        self.deactivate_elements(self.elements.length < 1, "Iceberg melted")
        self.deactivate_elements(self.elements.width < 1, "Iceberg melted")
        self.deactivate_elements(self.elements.sail < 1, "Iceberg melted")


    def roll_over(self, rho_water):
        """ Iceberg's stability criterium: Parameterization from Keghouche et al. 2009 with a correction from Wagner et al. 2017 """
        L = self.elements.length
        W = self.elements.width
        H = self.elements.draft + self.elements.sail
        alpha = rho_iceb / rho_water
        crit = np.sqrt(6 * alpha * (1 - alpha))
        W, L = np.min([L, W], axis=0), np.max([L, W], axis=0)
        mask = (W / H) < crit
        if any(mask):
            logger.info(f"Rolling over : {np.sum(mask)} icebergs ...")
            nL, nW, nH = (np.max([L[mask], H[mask]], axis=0), np.min([L[mask], H[mask]], axis=0), W[mask])
            L[mask], W[mask], H[mask] = nL, nW, nH
        depthib = H * alpha
        sailib = H - depthib
        self.elements.length = L
        self.elements.width = W
        self.elements.sail = sailib
        self.elements.draft = depthib

    def prepare_run(self):
        self.profiles_depth = self.elements_scheduled.draft.max()
        logger.info(f"Icebergs max draft is : {self.profiles_depth} meters")

    def update(self):
        """ Update positions and properties of particles """
        T = self.environment.sea_water_temperature
        S = self.environment.sea_water_salinity
        rho_water = PhysicsMethods.sea_water_density(T, S)
        self.roll_over(rho_water)
        if self.melting:
            self.melt()
        self.advect_iceberg(rho_water, self.add_stokes_drift, self.wave_rad, self.grounding, self.vertical_profile)