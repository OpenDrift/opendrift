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

from datetime import datetime, timedelta
import numpy as np
from scipy.interpolate import interp1d
import logging; logger = logging.getLogger(__name__)
from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements import LagrangianArray
from opendrift.models.physics_methods import verticaldiffusivity_Large1994, verticaldiffusivity_Sundby1983, gls_tke, skillscore_liu_weissberg
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED

# Defining the oil element properties
class Lagrangian3DArray(LagrangianArray):
    """Extending LagrangianArray for elements moving in 3 dimensions
    The Particle may be buoyant and/or subject to vertical mixing
    buoyant bahaviour is described by terminal velocity
    """

    variables = LagrangianArray.add_variables([
        ('wind_drift_factor', {'dtype': np.float32,
                               'units': '1',
            'description': 'Elements at surface are moved with this '
                'fraction of the vind vector, in addition to currents '
                'and Stokes drift',
                               'default': 0.02}),
        ('current_drift_factor', {'dtype': np.float32,
                                  'units': '1',
            'description': 'Elements are moved with this fraction of the '
                            'current vector, in addition to currents '
                            'and Stokes drift',
                               'default': 1}),
        ('terminal_velocity', {'dtype': np.float32,
                               'units': 'm/s',
            'description': 'Terminal rise/sinking velocity (buoyancy) '
                'in the ocean column',
                               'default': 0.})])


class OceanDrift(OpenDriftSimulation):
    """Open source buoyant particle trajectory model based on OpenDrift.

        Developed at MET Norway

        Generic module for particles that move in 3 dimensions
        and may be to vertical turbulent mixing
        with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, nutrients or sediments,
        Model may be subclassed for more specific behaviour.

    """

    ElementType = Lagrangian3DArray

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_height': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0},
        'ocean_vertical_diffusivity': {'fallback': 0,
             'profiles': True},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
        'sea_surface_wave_period_at_variance_spectral_density_maximum':
            {'fallback': 0},
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': {'fallback': 0},
        'sea_surface_swell_wave_to_direction': {'fallback': 0, 'important': False},
        'sea_surface_swell_wave_peak_period_from_variance_spectral_density': {'fallback': 0, 'important': False},
        'sea_surface_swell_wave_significant_height': {'fallback': 0, 'important': False},
        'sea_surface_wind_wave_to_direction': {'fallback': 0, 'important': False},
        'sea_surface_wind_wave_mean_period': {'fallback': 0, 'important': False},
        'sea_surface_wind_wave_significant_height': {'fallback': 0, 'important': False},
        'surface_downward_x_stress': {'fallback': 0},
        'surface_downward_y_stress': {'fallback': 0},
        'turbulent_kinetic_energy': {'fallback': 0},
        'turbulent_generic_length_scale': {'fallback': 0},
        'ocean_mixed_layer_thickness': {'fallback': 50},
        'sea_floor_depth_below_sea_level': {'fallback': 10000},
        'land_binary_mask': {'fallback': None},
        }


    def __init__(self, *args, **kwargs):

        if 'machine_learning_dict' in kwargs:
            logger.debug('Machine learning correction supplied.')
            mld = kwargs['machine_learning_dict']
            del kwargs['machine_learning_dict']
            import pickle
            from tensorflow import keras
            from tools.normalizations import generate_one_normalized_predictor, generate_one_denormalized_prediction
            self.generate_one_normalized_predictor = generate_one_normalized_predictor
            self.generate_one_denormalized_prediction = generate_one_denormalized_prediction
            self.trained_model = keras.models.load_model(mld['trained_model'])
            self.dict_normalization_params = pickle.load(open(mld['dict_normalization_params'], 'rb'))
            self.ml_predictors = [self.dict_normalization_params['predictors'][i]['content'] for i in self.dict_normalization_params['predictors']]
        else:
            logger.debug('No machine learning correction available.')

        if hasattr(self, 'required_profiles_z_range'):
            raise ValueError('self.required_profiles_z_range is obsolete, and replaced by config setting drift:profile_depth. Please update your model.')

        # Calling general constructor of parent class
        super(OceanDrift, self).__init__(*args, **kwargs)

        self._add_config({
            'drift:vertical_advection': {'type': 'bool', 'default': True, 'description':
                'Advect elements with vertical component of ocean current.',
                'level': CONFIG_LEVEL_BASIC},
            'drift:vertical_mixing': {'type': 'bool', 'default': False, 'level': CONFIG_LEVEL_BASIC,
                'description': 'Activate vertical mixing scheme with inner loop'},
            'vertical_mixing:timestep': {'type': 'float', 'min': 0.1, 'max': 3600, 'default': 60,
                'level': CONFIG_LEVEL_ADVANCED, 'units': 'seconds', 'description':
                'Time step used for inner loop of vertical mixing.'},
            'vertical_mixing:diffusivitymodel': {'type': 'enum', 'default': 'environment',
                'enum': ['environment', 'stepfunction', 'windspeed_Sundby1983',
                 'windspeed_Large1994', 'gls_tke','constant'], 'level': CONFIG_LEVEL_ADVANCED,
                 'units': 'seconds', 'description': 'Algorithm/source used for profile of vertical diffusivity. Environment means that diffusivity is aquired from readers or environment constants/fallback.'},
            'vertical_mixing:background_diffusivity': {'type': 'float', 'min': 0, 'max': 1, 'default': 1.2e-5,
                'level': CONFIG_LEVEL_ADVANCED, 'units': 'm2s-1', 'description':
                'Background diffusivity used below mixed layer for wind parameterisations.'},
            'vertical_mixing:TSprofiles': {'type': 'bool', 'default': False, 'level':
                CONFIG_LEVEL_ADVANCED,
                'description': 'Update T and S profiles within inner loop of vertical mixing. This takes more time, but may be slightly more accurate.'},
            'drift:wind_drift_depth': {'type': 'float', 'default': 0.1,
                'min': 0, 'max': 10, 'units': 'meters',
                'description': 'The direct wind drift (windage) is linearly decreasing from the surface value (wind_drift_factor) until 0 at this depth.',
                'level': CONFIG_LEVEL_ADVANCED},
            'drift:stokes_drift': {'type': 'bool', 'default': True,
                'description': 'Advection elements with Stokes drift (wave orbital motion).',
                'level': CONFIG_LEVEL_ADVANCED},
            'drift:stokes_drift_profile': {'type': 'enum', 'default': 'Phillips',
                                           'enum': ['monochromatic', 'exponential', 'Phillips', 'windsea_swell'],
                                           'description': 'Algorithm to calculate Stokes drift at depth from surface value',
                                           'level': CONFIG_LEVEL_ADVANCED},
            'drift:use_tabularised_stokes_drift': {'type': 'bool', 'default': False,
                'description': 'If True, Stokes drift is estimated from wind based on look-up-tables for given fetch (drift:tabularised_stokes_drift_fetch).',
                'level': CONFIG_LEVEL_ADVANCED},
            'drift:tabularised_stokes_drift_fetch': {'type': 'enum', 'enum': ['5000', '25000', '50000'], 'default': '25000',
                'level': CONFIG_LEVEL_ADVANCED, 'description':
                'The fetch length when using tabularised Stokes drift.'},
            'general:seafloor_action': {'type': 'enum', 'default': 'lift_to_seafloor',
                'enum': ['none', 'lift_to_seafloor', 'deactivate', 'previous'],
                'description': '"deactivate": elements are deactivated; "lift_to_seafloor": elements are lifted to seafloor level; "previous": elements are moved back to previous position; "none"; seafloor is ignored.',
                'level': CONFIG_LEVEL_ADVANCED},
            'drift:truncate_ocean_model_below_m': {'type': 'float', 'default': None,
                'min': 0, 'max': 10000, 'units': 'm',
                'description': 'Ocean model data are only read down to at most this depth, and extrapolated below. May be specified to read less data to improve performance.',
                'level': CONFIG_LEVEL_ADVANCED},
             'seed:z': {'type': 'float', 'default': 0,
                    'min': -10000, 'max': 0, 'units': 'm',
                'description': 'Depth below sea level where elements are released. This depth is neglected if seafloor seeding is set selected.',
                'level': CONFIG_LEVEL_ESSENTIAL},
            'seed:seafloor': {'type': 'bool', 'default': False,
                'description': 'Elements are seeded at seafloor, and seeding depth (z) is neglected.',
                'level': CONFIG_LEVEL_ESSENTIAL},
            })

        self._set_config_default('drift:max_speed', 2)

    def update(self):
        """Update positions and properties of elements."""

        # Simply move particles with ambient current
        self.advect_ocean_current()

        # Advect particles due to surface wind drag,
        # according to element property wind_drift_factor
        self.advect_wind()

        # Stokes drift
        self.stokes_drift()

        # Turbulent Mixing
        if self.get_config('drift:vertical_mixing') is True:
            self.update_terminal_velocity()
            self.vertical_mixing()
        else:  # Buoyancy
            self.vertical_buoyancy()

        # Vertical advection
        self.vertical_advection()

        # Optional machine learning correction
        self.machine_learning_correction()

    def simulate_trajectories(self, outfile, trajectories, number=1,
                              wind_drift_factors=None, current_drift_factors=None,
                              time_step=None, time_step_output=None,
                              simulation_duration=None, simulation_interval=None):

        import pandas as pd
        time_step_output = pd.Timedelta(time_step_output)
        simulation_interval = pd.Timedelta(simulation_interval)
        # Interpolate trajectories to output time step
        #trajectories = trajectories.traj.gridtime(time_step_output)
        # Find all seed combinations: position, time, wdf, cdf
        #  Loop også over trajektorier -> origin_marker
        step = int(simulation_interval / time_step_output)
        tind = np.arange(0, trajectories.sizes['time'], step)
        start_lons = trajectories.isel(time=tind).isel(trajectory=0).lon.values
        start_lats = trajectories.isel(time=tind).isel(trajectory=0).lat.values
        start_times = trajectories.time[tind]
        self.set_config('drift:max_age_seconds', simulation_duration.total_seconds())
        kwargs = {'wind_drift_factor': wind_drift_factors[0],
                  'current_drift_factor': current_drift_factors[0]}
        for (lo,la,ti) in zip(start_lons, start_lats, start_times):
            if np.isnan(lo):
                continue
            ti = pd.Timestamp(ti.values).to_pydatetime()
            self.seed_elements(lon=lo, lat=la, time=ti, number=number, **kwargs)
        print(self)
        self.run(outfile=outfile, end_time=pd.Timestamp(start_times[-1].values).to_pydatetime()+simulation_duration)
        # Simulate and save to file

    def wind_drift_factor_from_trajectory_lw(self, drifters, wind_drift_factors,
                                             simulation_length, simulation_interval):
        """Perform simulations and use skillscore to optimize wind_drift_factor

        drifters: list of dictionaries with numpy arrays of 'lon' and 'lat'
                    and list of datetimes
        wind_drift_factors: the wind_drift_factors to use for simulations/optimalizations
        """

        output = []  # List with one dictionary per trajectory
        for d in drifters:
            time = np.array(d['time'])
            ti = time[1::] - time[0:-1]
            if ti.min() == ti.max():
                trajectory_interval = ti[0]
            else:
                raise ValueError('Trajectory interval is not constant')
            if 'tiv' not in locals():
                tiv = trajectory_interval
                index_interval = simulation_interval.total_seconds()/trajectory_interval.total_seconds()
                if not index_interval.is_integer():
                    raise ValueError('Simulation interval must be a multiple of trajectory time step')
            else:
                if trajectory_interval != tiv:
                    raise ValueError('Trajectories do not have the same time step')
            last_index = (time[-1]-time[0]-simulation_length).total_seconds()/trajectory_interval.total_seconds()
            seed_indices = np.arange(0, last_index, index_interval).astype(int)
            do = {'lon': d['lon'][seed_indices],
                  'lat': d['lat'][seed_indices],
                  'time': np.array(time[seed_indices]),
                  'trajectory_start_index': seed_indices}
            output.append(do)

            if not 'last_seed_time' in locals():
                last_seed_time = do['time'][-1]
            else:
                last_seed_time = np.maximum(do['time'][-1], last_seed_time)
            # Seed for all starting positions in this trajectory
            for lo,la,ti in zip(do['lon'], do['lat'], do['time']):
                ow = np.ones(len(wind_drift_factors))
                self.seed_elements(lon=lo*ow, lat=la*ow, time=ti,
                                   wind_drift_factor=wind_drift_factors)

        self.set_config('drift:max_age_seconds', simulation_length.total_seconds()+1)
        self.run(end_time=last_seed_time+simulation_length)

        index_of_first, index_of_last = \
            self.index_of_activation_and_deactivation()
        times_model = np.array(self.get_time_array()[0])

        #import matplotlib.pyplot as plt
        i = 0  # element number in simulation
        for d, o in zip(drifters, output):
            o['segments'] = {}
            o['wind_drift_factor'] = np.ones(len(o['trajectory_start_index']))
            o['skillscore'] = np.ones(len(o['trajectory_start_index']))
            for segnum, tsi in enumerate(o['trajectory_start_index']):
                o['segments'][segnum] = {'skillscore': ow}
                lon_obs = d['lon'][tsi:tsi+int(index_interval)]
                lat_obs = d['lat'][tsi:tsi+int(index_interval)]
                #plt.plot(lon_obs, lat_obs, 'r')
                for wnum, wdf in enumerate(wind_drift_factors):
                    lon_model = self.history['lon'][i, index_of_first[i]:index_of_last[i]]
                    lat_model = self.history['lat'][i, index_of_first[i]:index_of_last[i]]
                    i = i + 1
                    ss = skillscore_liu_weissberg(lon_obs, lat_obs, lon_model, lat_model)
                    o['segments'][segnum]['skillscore'][wnum] = ss
                    #plt.plot(lon_model, lat_model, label='Skillscore: %s' % ss)

                print(o['segments'][segnum]['skillscore'], 'SS')
                o['wind_drift_factor'][segnum] = \
                    wind_drift_factors[np.argmax(o['segments'][segnum]['skillscore'])]
                o['skillscore'][segnum] = o['segments'][segnum]['skillscore'].max()
                #plt.legend()
                #plt.show()
                #stop

        #import xarray as xr
        #da = xr.DataArray(
        #    data=skillscore,
        #    dims=["wind_drift_factor", "time"],
        #    coords=dict(
        #        wind_drift_factor=wind_drift_factors,
        #        time=time
        #    ),
        #    attrs=dict(
        #        description="Liu Weissberg skillscore",
        #    ),
        #)

        return output

    def machine_learning_correction(self):
        if not hasattr(self, 'trained_model'):
            return  # No machine learning correction available

        # ML shall only be applied every 1 hour
        if not hasattr(self, 'ml_timesteps'):
            self.ml_timesteps = 3600/self.time_step.total_seconds()
        if self.steps_calculation % self.ml_timesteps != 0:
            logger.warning('ML not applied at time step %s' % self.steps_calculation)
            return

        # Only apply every 1 hour
        logger.warning('Machine learning')
        x_vel = np.ones(self.num_elements_active())*np.NaN
        y_vel = np.ones(self.num_elements_active())*np.NaN
        st = datetime.now()

        list_predictors_normalized = []
        list_predictors_denormalized = []

        # Prepare predictor dicts
        for n in range(self.num_elements_active()):
            dict_input = {s: getattr(self.environment, s)[n] for s in self.ml_predictors}
            predictor_normalized = self.generate_one_normalized_predictor(
                                            dict_input, self.dict_normalization_params)
            list_predictors_normalized.append(predictor_normalized)

        # Normalize
        all_predictors_normalized = np.squeeze(np.array(list_predictors_normalized))
        # perform a prediction from normalized predictor to normalized label
        logger.warning('Calculating ML correction...')
        predicted_normalized_residual_correction = self.trained_model.predict(all_predictors_normalized)
        # get the "native units" residual correction
        logger.warning('Denormalization...')
        for n in range(self.num_elements_active()):
            native_units_correction = self.generate_one_denormalized_prediction(
                predicted_normalized_residual_correction[n, :], self.dict_normalization_params)
            x_vel[n] = native_units_correction['residual_displacement_x']*1000/3600
            y_vel[n] = native_units_correction['residual_displacement_y']*1000/3600

        # Apply correction
        logger.warning('Applying ML correction: %s to %s m/s eastwards, %s to %s m/s northwards' %
                        (x_vel.min(), x_vel.max(), y_vel.min(), y_vel.max()))
        #logger.warning('%s particles, %s' % (self.num_elements_active(), datetime.now()-st))
        self.update_positions(x_vel, y_vel)

    def disable_vertical_motion(self):
        """Deactivate any vertical processes/advection"""
        conf = {
                'drift:vertical_advection': False,
                'drift:vertical_mixing': False}
        for co, va in conf.items():
            logger.info('Setting config: %s -> %s' % (co, va))
            self.set_config(co, va)

    def update_terminal_velocity(self, Tprofiles=None, Sprofiles=None,
                                 z_index=None):
        """Calculate terminal velocity due to bouyancy from own properties
        and environmental variables. Sub-modules should overload
        this method for particle-specific behaviour
        """
        pass

    def prepare_vertical_mixing(self):
        pass  # To be implemented by subclasses as needed

    def prepare_run(self):
        super(OceanDrift, self).prepare_run()

    def vertical_advection(self):
        """Move particles vertically according to vertical ocean current

            Vertical advection by ocean currents is normally small
            compared to termical velocity
        """
        if self.get_config('drift:vertical_advection') is False:
            logger.debug('Vertical advection deactivated')
            return

        in_ocean = np.where(self.elements.z<0)[0]
        if len(in_ocean) > 0:
            w = self.environment.upward_sea_water_velocity[in_ocean]
            self.elements.z[in_ocean] = np.minimum(0,
                self.elements.z[in_ocean] + self.elements.moving[in_ocean] * w * self.time_step.total_seconds())
        else:
            logger.debug('No vertical advection for elements at surface')

    def vertical_buoyancy(self):
        """Move particles vertically according to their buoyancy"""
        in_ocean = np.where(self.elements.z<0)[0]
        if len(in_ocean) > 0:
            self.elements.z[in_ocean] = np.minimum(0,
                self.elements.z[in_ocean] + self.elements.terminal_velocity[in_ocean] * self.time_step.total_seconds())

        # check for minimum height/maximum depth for each particle accouting also for
        # the sea surface height
        Zmin = -1.*(self.environment.sea_floor_depth_below_sea_level + self.environment.sea_surface_height)

        # Let particles stick to bottom
        bottom = np.where(self.elements.z < Zmin)
        if len(bottom[0]) > 0:
            logger.debug('%s elements reached seafloor, interacting with bottom' % len(bottom[0]))
            self.interact_with_seafloor()
            self.bottom_interaction(Zmin)

    def surface_stick(self):
        '''To be overloaded by subclasses, e.g. downward mixing of oil'''

        # keep particle just below the surface
        surface = np.where(self.elements.z >= 0)
        if len(surface[0]) > 0:
            self.elements.z[surface] = -0.01

    def bottom_interaction(self, Zmin=None):
        '''To be overloaded by subclasses, e.g. radionuclides in sediments'''
        pass

    def surface_wave_mixing(self, time_step_seconds):
        '''To be overloaded by subclasses, e.g. downward mixing of oil'''

        pass

    def get_diffusivity_profile(self, model, depths):
        wind, depth = np.meshgrid(self.wind_speed(), depths)
        MLD = self.environment.ocean_mixed_layer_thickness
        background_diffusivity = self.get_config('vertical_mixing:background_diffusivity')

        if model == 'windspeed_Large1994':
            return verticaldiffusivity_Large1994(wind, depth, MLD, background_diffusivity)
        elif model == 'windspeed_Sundby1983':
            return verticaldiffusivity_Sundby1983(wind, depth, MLD, background_diffusivity)
        elif model == 'gls_tke':
            if not hasattr(self, 'gls_parameters'):
                logger.info('Searching readers for GLS parameters...')
                for reader_name, reader in self.readers.items():
                    if hasattr(reader, 'gls_parameters'):
                        self.gls_parameters = reader.gls_parameters
                        logger.info('Found gls-parameters in ' + reader_name)
                        break  # Success
                if not hasattr(self, 'gls_parameters'):
                    logger.info('Did not find gls-parameters in any readers.')
                    self.gls_parameters = None
            windstress = np.sqrt(self.environment.surface_downward_x_stress**2 +
                                 self.environment.surface_downward_y_stress**2)
            return gls_tke(windstress, depth, self.sea_water_density(),
                           self.environment.turbulent_kinetic_energy,
                           self.environment.turbulent_generic_length_scale,
                           gls_parameters)

        else:
            raise ValueError('Unknown diffusivity model: ' + model)

    def vertical_mixing(self, store_depths=False):
        """Mix particles vertically according to eddy diffusivity and buoyancy

            Buoyancy is expressed as terminal velocity, which is the
            steady-state vertical velocity due to positive or negative
            buoyant behaviour. It is usually a function of particle density,
            diameter, and shape.

            Vertical particle displacemend du to turbulent mixing is
            calculated using a random walk scheme" (Visser et al. 1996)
        """

        if self.get_config('drift:vertical_mixing') is False:
            logger.debug('Turbulent mixing deactivated')
            return

        self.timer_start('main loop:updating elements:vertical mixing')

        dt_mix = self.get_config('vertical_mixing:timestep')

        # minimum height/maximum depth for each particle accouting also for
        # the sea surface height
        Zmin = -1.*(self.environment.sea_floor_depth_below_sea_level + self.environment.sea_surface_height)

        # Eventual model specific preparions
        self.prepare_vertical_mixing()

        # get profile of eddy diffusivity
        # get vertical eddy diffusivity from environment or specific model
        diffusivity_model = self.get_config('vertical_mixing:diffusivitymodel')
        diffusivity_fallback = self.get_config('environment:fallback:ocean_vertical_diffusivity')
        # Mixing levels to be used if analytical expression:
        mixing_z_analytical = -np.arange(0, self.environment.ocean_mixed_layer_thickness.max() + 2)
        if diffusivity_model == 'environment':
            if 'ocean_vertical_diffusivity' in self.environment_profiles and not (
                    self.environment_profiles['ocean_vertical_diffusivity'].min() ==
                    diffusivity_fallback and
                    self.environment_profiles['ocean_vertical_diffusivity'].max() ==
                    diffusivity_fallback):
                Kprofiles = self.environment_profiles[
                    'ocean_vertical_diffusivity']
                mixing_z = self.environment_profiles['z'].copy()
                logger.debug('Using diffusivity from ocean model')
            else:
                logger.debug('Using diffusivity from Large1994 since model diffusivities not available')
                mixing_z = mixing_z_analytical
                Kprofiles = self.get_diffusivity_profile('windspeed_Large1994', np.abs(mixing_z))
        elif diffusivity_model == 'constant':
            logger.debug('Using constant diffusivity specified by fallback_values[''ocean_vertical_diffusivity''] = %s m2.s-1' % (diffusivity_fallback))
            mixing_z = self.environment_profiles['z'].copy()
            Kprofiles = diffusivity_fallback*np.ones(
                    self.environment_profiles['ocean_vertical_diffusivity'].shape) # keep constant value for ocean_vertical_diffusivity
        else:
            logger.debug('Using functional expression for diffusivity')
            # Note: although analytical functions, z is discretised
            mixing_z = mixing_z_analytical
            Kprofiles = self.get_diffusivity_profile(diffusivity_model, np.abs(mixing_z))

        logger.debug('Diffusivities are in range %s to %s' %
                      (Kprofiles.min(), Kprofiles.max()))

        # get profiles of salinity and temperature
        # (to save interpolation time in the inner loop)
        if (self.get_config('vertical_mixing:TSprofiles') is True
            and 'sea_water_salinity' in self.required_variables):
            Sprofiles = self.environment_profiles['sea_water_salinity']
            Tprofiles = \
                self.environment_profiles['sea_water_temperature']
            if (self.get_config('environment:fallback:sea_water_salinity') is not None and
                Sprofiles.min() == Sprofiles.max()):
                logger.debug('Salinity and temperature are fallback'
                              'values, skipping TSprofile')
                Sprofiles = None
                Tprofiles = None
            else:
                logger.debug('Using TSprofiles for vertical mixing')
        else:
            logger.debug('TSprofiles deactivated for vertical mixing')
            Sprofiles = None
            Tprofiles = None

        # prepare vertical interpolation coordinates
        z_i = range(mixing_z.shape[0])
        if len(z_i) == 1:
            z_index = 0
        else:
            z_index = interp1d(-mixing_z, z_i, bounds_error=False,
                               fill_value=(0,len(z_i)-1))  # Extrapolation

        # Internal loop for fast time step of vertical mixing model.
        # Random walk needs faster time step than horizontal advection.
        logger.debug('Vertical mixing module:' +
            self.get_config('vertical_mixing:diffusivitymodel'))
        ntimes_mix = np.abs(int(self.time_step.total_seconds()/dt_mix))
        logger.debug('Turbulent diffusion with random walk '
                      'scheme using ' + str(ntimes_mix) +
                      ' fast time steps of dt=' + str(dt_mix) + 's')

        if store_depths is not False:
            depths = np.zeros((ntimes_mix, self.num_elements_active()))
            depths[0, :] = self.elements.z

        # Calculating dK/dz for all profiles before the loop
        gradK = -np.gradient(Kprofiles, mixing_z, axis=0)
        gradK[np.abs(gradK)<1e-10] = 0

        for i in range(0, ntimes_mix):
            #remember which particles belong to the exact surface
            surface = self.elements.z == 0

            # Update the terminal velocity of particles
            self.update_terminal_velocity(Tprofiles=Tprofiles, Sprofiles=Sprofiles, z_index=z_index)
            w = self.elements.terminal_velocity

            # Diffusivity and its gradient at z
            zi = np.round(z_index(-self.elements.z)).astype(np.uint16)
            Kz = Kprofiles[zi, range(Kprofiles.shape[1])]
            dKdz = gradK[zi, range(Kprofiles.shape[1])]

            # Visser et al. 1997 random walk mixing
            # requires an inner loop time step dt such that
            # dt << (d2K/dz2)^-1, e.g. typically dt << 15min
            #
            # NB: In the last term Kz is evaluated in zi, while
            # it should be evaluated in (self.elements.z - dKdz*dt_mix)
            # This is not expected have large impact on the result
            R = 2*np.random.random(self.num_elements_active()) - 1
            r = 1.0/3
            # New position  =  old position   - up_K_flux   + random walk
            self.elements.z = self.elements.z - self.elements.moving*(
                dKdz*dt_mix - R*np.sqrt((Kz*dt_mix*2/r)))

            # Reflect from surface
            reflect = np.where(self.elements.z >= 0)
            if len(reflect[0]) > 0:
                self.elements.z[reflect] = -self.elements.z[reflect]

            # Reflect elements going below seafloor
            # Zmin accounts for sea_surface_height
            bottom = np.where(np.logical_and(self.elements.z < Zmin, self.elements.moving == 1))
            if len(bottom[0]) > 0:
                logger.debug('%s elements penetrated seafloor, lifting up' % len(bottom[0]))
                self.elements.z[bottom] = 2*Zmin[bottom] - self.elements.z[bottom]

            # Advect due to buoyancy
            self.elements.z = self.elements.z + w*dt_mix*self.elements.moving

            # Put the particles that belonged to the surface slick
            # (if present) back to the surface
            self.elements.z[surface] = 0.

            # Formation of slick and wave mixing for surfaced particles
            # if implemented for this class
            self.surface_stick()
            self.surface_wave_mixing(dt_mix)

            # Let particles stick to bottom
            bottom = np.where(self.elements.z < Zmin)
            if len(bottom[0]) > 0:
                logger.debug('%s elements reached seafloor, interacting with bottom' % len(bottom[0]))
                self.interact_with_seafloor()
                self.bottom_interaction(Zmin)

            if store_depths is not False:
                depths[i, :] = self.elements.z

        self.timer_end('main loop:updating elements:vertical mixing')

        if store_depths is not False:
            return depths
        else:
            return None

    def animate_vertical_distribution(self, depths=None, maxdepth=None, bins=50, filename=None, subsamplingstep=1):
        """Function to animate vertical distribution of particles
            bins:            number of bins in the histogram
            maxdepth:        maximum depth
            subsamplingstep: speed-up the generation of the animation reducing the number of output frames
            fasterwriter:    speed-up the writing to outpute file
        """
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation

        start_time = datetime.now()

        #from timeit import default_timer as timer
        #start=timer()

        fig, (axk, axn) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3]})
        if depths is not None:  # Debug mode, output from one cycle has been provided
            z = depths
            time_step = self.get_config('vertical_mixing:timestep')
            times = [self.time + i*timedelta(seconds=time_step) for i in range(z.shape[0])]
        else:
            z = self.get_property('z')[0]
            z = np.ma.filled(z, np.nan)
            K = self.get_property('ocean_vertical_diffusivity')[0]
            time_step = self.time_step.total_seconds()
            times = self.get_time_array()[0]
        if maxdepth is None:
            maxdepth = np.nanmin(z)
        if maxdepth > 0:
            maxdepth = -maxdepth  # negative z

        if depths is not None:
            axk.plot(np.nanmean(self.environment_profiles['ocean_vertical_diffusivity'], 1), self.environment_profiles['z'])
            xmax = self.environment_profiles['ocean_vertical_diffusivity'].max()
        else:
            axk.plot(K, z, 'k.')
            xmax = np.nanmax(K)
        axk.set_ylim([maxdepth, 0])
        axk.set_xlim([0, xmax*1.1])
        axk.set_ylabel('Depth [m]')
        axk.set_xlabel('Vertical diffusivity [$m^2/s$]')

        subsamplingstep=max(1,subsamplingstep)
        times=times[0:-1:subsamplingstep]
        z=z[0:-1:subsamplingstep,:]

        hist_series = np.zeros((bins, len(times)))
        bin_series = np.zeros((bins+1, len(times)))
        for i in range(len(times)):
            hist_series[:,i], bin_series[:,i] = np.histogram(z[i,:][np.isfinite(z[i,:])], bins=bins, range=(maxdepth,z[np.isfinite(z)].max()))
        maxnum = hist_series.max()

        axn.clear()
        bc=axn.barh(bin_series[0:-1,i], hist_series[:,i], height=-maxdepth/bins, align='edge')
        axn.set_ylim([maxdepth, 0])
        axn.set_xlim([0, maxnum])
        title=axn.set_title('%s UTC' % times[i])
        axn.set_xlabel('Number of particles')

        def update_histogram(i):
            for rect, y in zip(bc, hist_series[:,i]):
                rect.set_width(y)
            title.set_text('%s UTC' % times[i])

        self.__save_or_plot_animation__(plt.gcf(),
                                        update_histogram,
                                        filename,
                                        frames = len(times),
                                        fps = 10,
                                        interval=50,
                                        blit=False)

        logger.info('Time to make animation: %s' %
                    (datetime.now() - start_time))




    def plot_vertical_distribution(self, maxdepth=None, bins=None, maxnum=None):
        """Function to plot vertical distribution of particles

            maxdepth: maximum depth considered for the profile
            bins:     number of bins between surface and maxdepth
            maxnum:   range of bars in histogram is [0,maxnum]
        """
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Slider

        fig = plt.figure()
        mainplot = fig.add_axes([.15, .3, .8, .5])
        sliderax = fig.add_axes([.15, .08, .75, .05])
        tslider = Slider(sliderax, 'Timestep', 0, self.steps_output-1,
                         valinit=self.steps_output-1, valfmt='%0.0f')
        try:
            dz = self.get_config('vertical_mixing:verticalresolution')
        except:
            dz = 1.
        maxrange = -100

        # overwrite default values if input arguments are provided
        if maxdepth is not None:
            maxrange = maxdepth
        if maxrange > 0:
            maxrange = -maxrange  # negative z

        if bins is not None:
            dz = -maxrange/bins

        # using get_property instead of history to exclude elements thare are not yet seeded
        z = self.get_property('z')[0]
        z = np.ma.filled(z, np.nan)

        if maxnum is None:
            # Precalculatig histograms to find maxnum
            hist_series = np.zeros((int(-maxrange/dz), self.steps_output-1))
            bin_series = np.zeros((int(-maxrange/dz)+1, self.steps_output-1))
            for i in range(self.steps_output-1):
                hist_series[:,i], bin_series[:,i] = np.histogram(z[i,:][np.isfinite(z[i,:])], bins=int(-maxrange/dz), range=[maxrange, 0])
            maxnum = hist_series.max()

        def update(val):
            tindex = int(tslider.val)
            mainplot.cla()
            mainplot.grid()
            mainplot.hist(z[tindex, :], bins=int(-maxrange/dz),
                          range=[maxrange, 0], orientation='horizontal')
            mainplot.set_ylim([maxrange, 0])
            mainplot.set_xlim([0, maxnum])
            mainplot.set_xlabel('number of particles')
            mainplot.set_ylabel('depth [m]')
            x_wind = self.history['x_wind'].T[tindex, :]
            y_wind = self.history['y_wind'].T[tindex, :]
            windspeed = np.mean(np.sqrt(x_wind**2 + y_wind**2))
            mainplot.set_title(str(self.get_time_array()[0][tindex]) +
                               #'   Percent at surface: %.1f %' % percent_at_surface)
                               '   Mean windspeed: %.1f m/s' % windspeed)
            fig.canvas.draw_idle()

        update(0)  # Plot initial distribution
        tslider.on_changed(update)
        plt.show()

        # returning objects prevents unresponsiveness when moving the Slider
        # https://github.com/matplotlib/matplotlib/issues/3105#issuecomment-44855888
        return fig,mainplot,sliderax,tslider

    def plotter_vertical_distribution_time(self, ax=None, mask=None,
            dz=1., maxrange=-100, bins=None, step=1):
        """Function to plot vertical distribution of particles.

        Use mask to plot any selection of particles.
        """
        from pylab import axes, draw
        from matplotlib import dates, pyplot

        if ax is None:
            fig = pyplot.figure()
            ax = fig.gca()
            show = True
        else:
            show = False

        if mask is None: # create a mask that is True for all particles
            mask = self.history['z'].T[0] == self.history['z'].T[0]

        if bins is None:
            bins=int(-maxrange/dz)

        ax.hist(self.history['z'].T[step,mask], bins=bins,
                range=[maxrange, 0], orientation='horizontal')
        ax.set_ylim([maxrange, 0])
        ax.grid()
        #ax.set_xlim([0, mask.sum()*.15])
        ax.set_xlabel('Number of particles')
        ax.set_ylabel('Depth [m]')
        x_wind = self.history['x_wind'].T[step, :]
        y_wind = self.history['x_wind'].T[step, :]
        windspeed = np.mean(np.sqrt(x_wind**2 + y_wind**2))
        ax.set_title(str(self.get_time_array()[0][step]) +
                     '   Mean windspeed: %.1f m/s' % windspeed)
        if show is True:
            pyplot.show()
