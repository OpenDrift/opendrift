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
# Copyright 2015, Knut-Frode Dagestad, MET Norway

# This ShipDrift model is based on
# Soergaard, E and T Vada (1998);
# Observations and modelling of drifting ships. Report DnV 96-2011
# Reprogrammed in Python for OpenDrift by Knut-Frode Dagestad Dec 2016

import os
import numpy as np
import scipy
import logging; logger = logging.getLogger(__name__)

from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements import LagrangianArray
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED


class ShipObject(LagrangianArray):
    """Extending LagrangianArray with variables relevant for leeway objects.

    """
    variables = LagrangianArray.add_variables([
        ('orientation', {'dtype': np.uint8,
                         'units': '1',
                         'default': 1}),
        ('length', {'dtype': np.float32,
                         'units': 'm',
                         'min': 1,
                         'max': 500,
            'description': 'Length of ship',
            'level': CONFIG_LEVEL_ESSENTIAL,
                         'default': 80}),
        ('height', {'dtype': np.float32,
                         'units': 'm',
                         'min': 1,
                         'max': 100,
            'description': 'Total height of ship (above and below waterline)',
            'level': CONFIG_LEVEL_ESSENTIAL,
                         'default': 8}),
        ('draft', {'dtype': np.float32,  # wet part of ship [m]
                         'units': 'm',
                         'min': 1,
                         'max': 30,
            'description': 'Draft of ship (depth below water)',
            'level': CONFIG_LEVEL_ESSENTIAL,
                         'default': 4.0}),
        ('beam', {'dtype': np.float32,  # width of ship
                         'min': 1,
                         'max': 70,
                         'units': 'm',
            'description': 'Beam (width) of ship',
            'level': CONFIG_LEVEL_ESSENTIAL,
                         'default': 10}),
        ('wind_drag_coeff', {'dtype': np.float32,  # Cf
                             'units': '1',
                             'default': 1}),
        ('water_drag_coeff', {'dtype': np.float32,  # Cd
                              'units': '1',
                              'default': 1}),
        ('jibeProbability', {'dtype': np.float32,
                             'units': '1/h',
                             'default': 0.04}),
       ])


class ShipDrift(OpenDriftSimulation):
    """Demonstration trajectory model based on OpenDrift framework.

    Simply advects a particle (passive tracer with
    no properties except for position) with the ambient wind.
    """

    ElementType = ShipObject

    required_variables = {
        'x_wind': {'fallback': None},
        'y_wind': {'fallback': None},
        'land_binary_mask': {'fallback': None},
        'x_sea_water_velocity': {'fallback': None},
        'y_sea_water_velocity': {'fallback': None},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': {'fallback': 0}
        }

    winwav_angle = 20  # Angular offset in degrees

    def __init__(self, *args, **kwargs):

        # Read ship properties
        d = os.path.dirname(os.path.realpath(__file__))
        w = open(d + '/wforce.dat', 'r')
        self.wforce = {}
        w.readline()
        nbeam = int(w.readline().split()[0])
        self.wforce['nbeam'] = nbeam
        self.wforce['BL'] = np.array(w.readline().split()[0:nbeam],
                                     dtype=float)
        ndraft = int(w.readline().split()[0])
        self.wforce['ndraft'] = ndraft
        self.wforce['DL'] = np.array(w.readline().split()[0:ndraft],
                                     dtype=float)
        nomega = int(w.readline().split()[0])
        self.wforce['nomega'] = nomega
        self.wforce['omega'] = np.zeros((nomega))
        self.wforce['F'] = np.zeros((nomega, nbeam, ndraft))
        self.wforce['D'] = np.zeros((nomega, nbeam, ndraft))
        for o in range(nomega):
            self.wforce['omega'][o] = float(w.readline().split()[0])
            for i in range(ndraft):
                l = w.readline()
                self.wforce['F'][o, i, :] = l.split()[0:nbeam]
            for i in range(ndraft):
                l = w.readline()
                self.wforce['D'][o, i, :] = l.split()[0:nbeam]

        w.close()
        wi_omega, wi_BL, wi_DL = \
            np.meshgrid(self.wforce['omega'],
                        self.wforce['BL'], self.wforce['DL'], indexing='ij')
        self.wforce_interpolator_F = scipy.interpolate.LinearNDInterpolator(
            (wi_omega.ravel(), wi_BL.ravel(), wi_DL.ravel()),
             self.wforce['F'].ravel())
        self.wforce_interpolator_D = scipy.interpolate.LinearNDInterpolator(
            (wi_omega.ravel(), wi_BL.ravel(), wi_DL.ravel()),
             self.wforce['D'].ravel())

        super(ShipDrift, self).__init__(*args, **kwargs)

        self._add_config({'seed:orientation': {'type': 'enum',
            'enum':['left', 'right', 'random'], 'default': 'random',
            'level': CONFIG_LEVEL_ESSENTIAL,
            'description': 'If ships are oriented to the left or right of the downwind direction,'
                'or whether this is unknown. Left/right means that wind will hit ship from backboard/steerboard'}})

        # Since the ShipDrift model is deterministic for given ship size,
        # (in contrast to the Leeway model), we use a default diffusivity
        # to yield some variability.
        self._set_config_default('drift:horizontal_diffusivity', 100)

        self._set_config_default('drift:max_speed', 2)

    def seed_elements(self, *args, **kwargs):

        if 'number' in kwargs:
            num = kwargs['number']
        else:
            num = self.get_config('seed:number')
        for var in ['length', 'height', 'draft', 'beam']:
            if var not in kwargs:
                kwargs[var] = self.get_config('seed:' + var)
            kwargs[var] = np.atleast_1d(kwargs[var])
            if len(kwargs[var] == 1):
                kwargs[var] = kwargs[var]*np.ones(num)

        # Check that beam and height vs length are within expected range
        dl = kwargs['draft']/kwargs['length']
        if dl.min() < 0.025 or dl.max() > 0.07:
            logger.warning('Ratio of draft to length should be in range '
                                '0.025 to 0.07, given range is %s-%s. '
                                'Using border value.' %
                                (dl.min(), dl.max()))
            dl = np.clip(dl, 0.025, 0.07)
        bl = kwargs['beam']/kwargs['length']
        if bl.min() < 0.12 or bl.max() > 0.18:
            logger.warning('Ratio of beam to length should be in range '
                                '0.12 to 0.18, given range is %s-%s. '
                                'Using border value.' %
                                (bl.min(), bl.max()))

        # Calculate drag coefficients based on given ship dimensions
        # Wind drag coefficient
        exposed = kwargs['height'] - kwargs['draft']
        Cf = np.zeros(num)
        Cf[exposed>37.2] = 1.4
        Cf[exposed<=37.2] = 1.045 + 0.016*(exposed[exposed<=37.2] - 15.)
        Cf[exposed<=15] = 0.700 + 0.023*exposed[exposed<=15]
        kwargs['wind_drag_coeff'] = Cf

        # Water drag coefficient
        beta = 2.0*dl
        Cd = np.zeros(num)
        Cd[beta>.12] = 1.27
        Cd[beta<=.12] =  1.32 + (1.27-1.32)/0.02 * (beta[beta<=.12]-0.10)
        Cd[beta<=.10] =  1.38 + (1.32-1.38)/0.02 * (beta[beta<=.10]-0.08)
        Cd[beta<=.08] =  1.44 + (1.38-1.44)/0.02 * (beta[beta<=.08]-0.06)
        Cd[beta<=.06] =  1.50 + (1.44-1.50)/0.01 * (beta[beta<=.06]-0.05)
        kwargs['water_drag_coeff'] = Cd

        if 'orientation' not in kwargs:
            oc = self.get_config('seed:orientation')
            if oc == 'left':
                kwargs['orientation'] = np.ones(num)*0
            elif oc == 'right':
                kwargs['orientation'] = np.ones(num)*1
            else:
                kwargs['orientation'] = np.r_[:num] % 2  # Random 0 or 1

        # Calling general constructor with calculated values
        super(ShipDrift, self).seed_elements(*args, **kwargs)

    def update(self):

        Tm = self.wave_period()
        Hs = self.significant_wave_height()
        dl = self.elements.draft/self.elements.length
        bl = self.elements.beam/self.elements.length
        # Clip to allowed range
        bl = np.clip(bl, 0.12, 0.18)
        dl = np.clip(dl, 0.025, 0.07)
        # Additional clipping to avoid NaN from interpolator
        bl = np.clip(bl, 0.121, 0.179)
        dl = np.clip(dl, 0.0251, 0.069)

        # Simply move particles with ambient current
        self.update_positions(self.environment.x_sea_water_velocity,
                              self.environment.y_sea_water_velocity)

        # Wind force
        rho_air = 1.25  # to be checked
        area_dry = self.elements.length*(self.elements.height -
                                         self.elements.draft)
        area_wet = self.elements.length*self.elements.draft
        # Wind force
        F_wind = (0.5*rho_air*self.elements.wind_drag_coeff*
                  area_dry*np.power(self.wind_speed(), 2))
        # Decompose wind
        F_wind_x = F_wind*self.environment.x_wind/self.wind_speed()
        F_wind_y = F_wind*self.environment.y_wind/self.wind_speed()
        F_wind_x[self.wind_speed()==0] = 0
        F_wind_y[self.wind_speed()==0] = 0

        # Wave force
        rho_water = 1025
        NSPEC = 100
        ommin2 = 2.25
        ommin3 = 7.0
        ommax = 12.0
        dom = (ommax-ommin2)/(NSPEC-1)
        F_wave = np.zeros(self.num_elements_active())  # For each ship
        beta2 = 0
        scale1 = np.sqrt(9.81/self.elements.length)

        tmp = np.power(2.0*np.pi/Tm, 4)
        d = tmp*Hs*Hs/(4*np.pi)
        b = tmp/np.pi
        s = np.zeros((NSPEC, self.num_elements_active()))
        for i in range(NSPEC):
            omi = ommin2 + i*dom
            omi = omi*scale1
            s[i, :] = d * np.exp(-b/np.power(omi, 4)) / np.power(omi, 5)  # m2s

        f2 = 0.0
        d2 = 0.0
        for i in range(NSPEC):
            omi = ommin2 + i*dom
            f1 = f2
            d1 = d2
            #print omi - ommin3, 'should be positive'
            #print omi
            if (omi < ommin3):
                f2 = self.wforce_interpolator_F(omi, bl, dl)
                d2 = self.wforce_interpolator_D(omi, bl, dl)
            else:
                # Interval 3
                f2 = 0.5
                d2 = 4.0*omi*f2

            F_wave = F_wave + 0.5*(f1+f2)*dom*scale1*np.power(s[i,:], 2)
            beta2 = beta2 + 0.5*(d1+d2)*dom*scale1*np.power(s[i,:], 2)

        F_wave = F_wave*rho_water*9.81*self.elements.length
        beta2 = beta2*rho_water*np.sqrt(9.81*self.elements.length)

        # Add calculated wave and wind drift
        longperiod = Tm > 8.55
        F_wave_b = F_wave.copy()
        F_wave[longperiod] = F_wave[longperiod]*.66
        beta2[longperiod] = beta2[longperiod]*.60
        medperiod = ((Tm >= 5.7) & (Tm <= 8.55))
        F_wave[medperiod] = F_wave[medperiod]*(1.0 - 0.34*(Tm[medperiod]-5.7)/2.85)
        beta2[medperiod] = beta2[medperiod]*(1.0 - 0.4*(Tm[medperiod]-5.7)/2.85)

        # Form drag (water resistance)
        beta1 = 0.5*rho_water*self.elements.water_drag_coeff*area_wet

        # Wave direction is taken as wind direction plus offset +/- 20 degrees
        # Using minus, since angles are here defined counter-clockwise from x-axis
        offset = -self.winwav_angle*2*(self.elements.orientation - 0.5)
        if (self.environment.sea_surface_wave_stokes_drift_x_velocity.max() == 0 and
            self.environment.sea_surface_wave_stokes_drift_y_velocity.max() == 0):
                logger.info('Using wind direction as wave direction')
                wave_dir = np.radians(offset) + np.arctan2(self.environment.y_wind,
                                                           self.environment.x_wind)
        else:
            logger.info('Using Stokes drift direction as wave direction')
            wave_dir = np.radians(offset) + np.arctan2(
                self.environment.sea_surface_wave_stokes_drift_y_velocity,
                self.environment.sea_surface_wave_stokes_drift_x_velocity)
        F_wave_x = F_wave*np.cos(wave_dir)
        F_wave_y = F_wave*np.sin(wave_dir)
        F_total = np.sqrt(np.power(F_wind_x + F_wave_x, 2) +
                          np.power(F_wind_y + F_wave_y, 2))

        # Iterate 4 times in order to estimate the effect of
        # wave damping and form drag
        uw_tot = 0
        uw_dir = 0
        for i in range(4):
            # Calc wave damping force in x and y (Ref (1), eq. 6.6.4)
            f2x = beta2*uw_tot*np.cos(wave_dir)
            f2y = beta2*uw_tot*np.sin(wave_dir)
            # Calc new drift direction, including effect of
            # wave damping (Ref (1), eq. 6.6.1)
            uw_dir = np.arctan2(F_wind_y+F_wave_y-f2y, F_wind_x+F_wave_x-f2x)
            # Ref (1), eq. 6.6.3
            bet2c = beta2*np.cos((wave_dir-uw_dir))
            # Ref (1), eq. 6.4.3
            uw_tot = -bet2c/(2.*beta1) + np.sqrt(bet2c*bet2c +
                                                 4.*beta1*F_total)/(2.*beta1)

        # Finally advect according to wind-wave forces
        velocity_u = uw_tot*np.cos(uw_dir)
        velocity_v = uw_tot*np.sin(uw_dir)
        self.update_positions(velocity_u, velocity_v)

        # Stranding
        self.deactivate_elements(self.environment.land_binary_mask == 1,
                                 reason='ship stranded')

    def import_ascii_format(self, filename):

        with open(filename) as f:
            lines = f.readlines()
        self.time_step_output, self.time_step = lines[16].split()
        num_elements = int(lines[25].split()[0])
        if num_elements != 1:
            raise ValueError('Import presently only supports single ship')
        num_timesteps = (len(lines)-30.)/14.
        num_timesteps = int(num_timesteps)

        # Initialise history array
        from datetime import datetime, timedelta
        history_dtype_fields = [
            (name, self.ElementType.variables[name]['dtype'])
            for name in self.ElementType.variables]
        # Add environment variables
        self.history_metadata = self.ElementType.variables.copy()
        for env_var in self.required_variables:
            history_dtype_fields.append((env_var, np.dtype('float32')))
            self.history_metadata[env_var] = {}
        history_dtype = np.dtype(history_dtype_fields)

        self.history = np.ma.array(
	        np.zeros([num_elements, num_timesteps]),
            dtype=history_dtype, mask=[True])

        self.steps_output = num_timesteps
        self.steps = num_timesteps
        self.start_time = datetime.now()
        self.end_time = datetime.now() + timedelta(hours=1)
        self.time = self.end_time

        # Read time steps from file
        for i in range(num_timesteps):
            line = lines[30 + 14*i + 13]
            l = line.split()
            lon = float(l[2])
            lat = float(l[3])
            self.history['lon'][0, i] = lon
            self.history['lat'][0, i] = lat
