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
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway

import numpy as np
import logging
from opendrift.models.opendrift3D import OpenDrift3DSimulation
from opendrift.models.oceandrift import OceanDrift
# from opendrift.elements.passivetracer import PassiveTracer
from opendrift.elements.buoyanttracer import BuoyantTracer
# from opendrift.elements.sedimenttracer import SedimentTracer


class SedimentDrift3D(OpenDrift3DSimulation): # based on OpenDrift3DSimulation base class
    """Trajectory model based on the OpenDrift framework using the OpenDrift3DSimulation baseclass

    Sediment 3D motion 
    Propagation with horizontal and vertical ocean currents, horizontal and 
    vertical diffusions (additional wind drag inherited from base class if needed).
    Suitable for sediment tracers, e.g. for tracking sediment particles.
    Adapted from OpenDrift3DSimulation by Simon Weppe - MetOcean Solutions.

    """
    ElementType = BuoyantTracer # simply use BuoyantTracer for now - will eventually move to SedimentTracer
    required_variables = [
        'x_sea_water_velocity',
        'y_sea_water_velocity',
        'x_wind', 'y_wind',
        'upward_sea_water_velocity',
        'ocean_vertical_diffusivity',
        'ocean_horizontal_diffusivity',
        'sea_surface_wave_significant_height',
        'sea_surface_wave_stokes_drift_x_velocity',
        'sea_surface_wave_stokes_drift_y_velocity',
        'sea_surface_wave_period_at_variance_spectral_density_maximum',
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
        'sea_floor_depth_below_sea_level'
        ]

    required_variables.append('land_binary_mask')

    required_profiles = ['ocean_vertical_diffusivity']
    # The depth range (in m) which profiles shall cover
    required_profiles_z_range = [-120, 0]

    fallback_values = {
        'x_sea_water_velocity': 0,
        'y_sea_water_velocity': 0,
        'sea_surface_wave_significant_height': 0,
        'sea_surface_wave_stokes_drift_x_velocity': 0,
        'sea_surface_wave_stokes_drift_y_velocity': 0,
        'sea_surface_wave_period_at_variance_spectral_density_maximum': 0,
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 0,
        'x_wind': 0,
        'y_wind': 0,
        'upward_sea_water_velocity': 0.0,
        'ocean_vertical_diffusivity': 0.0, # in m2/s
        'ocean_horizontal_diffusivity' : 0.0, # in m2/s2
        'sea_floor_depth_below_sea_level': 10000
        }
    
    # Adding some specs - inspired from basereader.py
    #
    # Default plotting colors of trajectory endpoints

    status_colors_default = {'initial': 'green',
                             'active': 'blue',
                             'missing_data': 'gray',
                             'settled': 'red'}
    # adding processes switch to contol resuspension
    configspec_sedimentdrift3d = '''
        [processes]
            resuspension = boolean(default=False)
    '''

    def __init__(self, *args, **kwargs):
        # update configstring
        self._add_configstring(self.configspec_sedimentdrift3d)

        # Calling general constructor of parent class
        super(SedimentDrift3D, self).__init__(*args, **kwargs)

        # Turbulent mixing switched off by default
        self.set_config('processes:turbulentmixing', False)
        # resuspension switched off by default
        self.set_config('processes:resuspension', False)

        self.max_speed = 1.5

    def get_environment(self, variables, time, lon, lat, z, profiles):
        '''Retrieve environmental variables at requested positions.
        This version was initially copied from basemodel.py, and modified to 
        allow possible addition of (x_sea_water_velocity,y_sea_water_velocity) 
        pairs, from several readers. 

        The intial version will just use a "break" after the variable group 
        (x_sea_water_velocity,y_sea_water_velocity) has been found in a reader
        (i.e. not adding current components)
        
        This is useful when tidal and residual currents come from different 
        datasets/files.


        Updates:
            Buffer (raw data blocks) for each reader stored for performance:
                [readers].var_block_before (last before requested time)
                [readers].var_block_after (first after requested time)
                    - lists of one ReaderBlock per variable group:
                        - time, x, y, [vars]
        Returns:
            environment: recarray with variables as named attributes,
                         interpolated to requested positions/time.

        '''
        add_water_velocity = False

        self.timer_start('main loop:readers')
        # Initialise ndarray to hold environment variables
        dtype = [(var, np.float32) for var in variables]
        env = np.ma.array(np.zeros(len(lon)), dtype=dtype)

        # For each variable/reader group:
        variable_groups, reader_groups, missing_variables = \
            self.get_reader_groups(variables)

        # Here, the function get_reader_groups
        # will take into account dataset priority
        for variable in variables:  # Fill with fallback value if no reader
            if (self.fallback_values is not None
                    and variable in self.fallback_values):
                env[variable] = np.ma.ones(env[variable].shape)\
                    * self.fallback_values[variable]

        for i, variable_group in enumerate(variable_groups):
            logging.debug('----------------------------------------')
            logging.debug('Variable group %s' % (str(variable_group)))
            logging.debug('----------------------------------------')
            reader_group = reader_groups[i]
            missing_indices = np.array(range(len(lon)))
            
            # -------------------------------------------------------------------------------
            # allow addition of several ['x_sea_water_velocity','y_sea_water_velocity'] pairs
            # add_water_velocity = False by default
            if 'x_sea_water_velocity' in variable_group and len(reader_group)>1:
            # i.e. the "current" variables interpolated are ['x_sea_water_velocity','y_sea_water_velocity']
            # and there are several readers for these : allow addition
                add_water_velocity = True 
            #--------------------------------------------------------------------------------

            # For each reader:
            for reader_name in reader_group:

                logging.debug('Calling reader ' + reader_name)
                logging.debug('----------------------------------------')
                self.timer_start('main loop:readers:' +
                                 reader_name.replace(':', '<colon>'))
                reader = self.readers[reader_name]
                if not reader.covers_time(time):
                    logging.debug('\tOutside time coverage of reader.')
                    continue
                # Fetch given variables at given positions from current reader
                try:
                    logging.debug('Data needed for %i elements' %
                                  len(missing_indices))
                    # Check if vertical profiles are requested from reader
                    if profiles is not None:
                        profiles_from_reader = list(
                            set(variable_group) & set(profiles))
                        if profiles_from_reader == []:
                            profiles_from_reader = None
                    else:
                        profiles_from_reader = None

                    if 'env_tmp' in locals() and add_water_velocity and 'x_sea_water_velocity' in reader.variables: 
                        # adding a special case when we add  ['x_sea_water_velocity','y_sea_water_velocity']
                        # this happens only if :
                        # env_tmp exists, which means a first interpoloation round happened
                        # there are more than 1 "reader" for the velocities, and
                        # we do want to use the velocity from that reader (possibly specified by variables_to_use when calling the reader)
                        #
                        logging.debug('Several readers available for %s : ' % (variable_group) )
                        logging.debug('%s' % (reader_group) )
                        logging.debug('Adding currents from reader : ' +  reader_name )
                        logging.debug('to existing [x_sea_water_velocity,y_sea_water_velocity]')
                        
                        # need to reset missing_indices to interpolate to all points
                        missing_indices = np.array(range(len(lon)))

                        # get (u,v) currents from that reader
                        env_tmp1, env_profiles_tmp1 = \
                            reader.get_variables_interpolated(
                                variable_group, profiles_from_reader,
                                self.required_profiles_z_range, time,
                                lon[missing_indices], lat[missing_indices],
                                z[missing_indices], self.use_block, self.proj)
                        
                        # add to previous currents
                        env_tmp['x_sea_water_velocity'] += env_tmp1['x_sea_water_velocity']
                        env_tmp['y_sea_water_velocity'] += env_tmp1['y_sea_water_velocity']

                        # add to previous current profile, if applicable
                        if env_profiles_tmp is not None : 
                            env_profiles_tmp['x_sea_water_velocity'] += env_profiles_tmp1['x_sea_water_velocity']
                            env_profiles_tmp['y_sea_water_velocity'] += env_profiles_tmp1['y_sea_water_velocity']

                    else :
                    # standard case - for all other variables groups or when there are only one reader for 
                    # ['x_sea_water_velocity','y_sea_water_velocity']
                    # Fetch given variables at given positions from current reader
                        env_tmp, env_profiles_tmp = \
                            reader.get_variables_interpolated(
                                variable_group, profiles_from_reader,
                                self.required_profiles_z_range, time,
                                lon[missing_indices], lat[missing_indices],
                                z[missing_indices], self.use_block, self.proj)


                except Exception as e:
                    logging.info('========================')
                    logging.info('Exception:')
                    logging.info(e)
                    logging.debug(traceback.format_exc())
                    logging.info('========================')
                    self.timer_end('main loop:readers:' +
                                   reader_name.replace(':', '<colon>'))
                    continue

                # Copy retrieved variables to env array, and mask nan-values
                for var in variable_group:
                    env[var][missing_indices] = np.ma.masked_invalid(
                        env_tmp[var][0:len(missing_indices)]).astype('float32')
                    if profiles_from_reader is not None and var in profiles_from_reader:
                        if 'env_profiles' not in locals():
                            env_profiles = env_profiles_tmp
                        # TODO: fix to be checked
                        if var in env_profiles and var in env_profiles_tmp:
                            # If one profile has fewer vertical layers than
                            # the other, we use only the overlapping part
                            if len(env_profiles['z']) != len(
                                env_profiles_tmp['z']):
                                logging.debug('Warning: different number of '
                                    ' vertical layers: %s and %s' % (
                                        len(env_profiles['z']),
                                        len( env_profiles_tmp['z'])))
                            z_ind = np.arange(np.minimum(
                                len(env_profiles['z'])-1,
                                len(env_profiles_tmp['z'])-1))
                            # len(missing_indices) since 2 points might have been added and not removed
                            env_profiles[var][np.ix_(z_ind, missing_indices)] = \
                                np.ma.masked_invalid(env_profiles_tmp[var][z_ind,0:len(missing_indices)]).astype('float32')
                
                # Detect elements with missing data, for present reader group
                if hasattr(env_tmp[variable_group[0]], 'mask'):
                    try:
                        del combined_mask
                    except:
                        pass
                    for var in variable_group:
                        tmp_var = np.ma.masked_invalid(env_tmp[var])
                        # Changed 13 Oct 2016, but uncertain of effect
                        # TODO: to be checked
                        #tmp_var = env_tmp[var]
                        if 'combined_mask' not in locals():
                            combined_mask = np.ma.getmask(tmp_var)
                        else:
                            combined_mask = \
                                np.ma.mask_or(combined_mask,
                                              np.ma.getmask(tmp_var),
                                              shrink=False)
                    try:
                        if len(missing_indices) != len(combined_mask):
                            # TODO: mask mismatch due to 2 added points
                            raise ValueError('Mismatch of masks')
                        missing_indices = missing_indices[combined_mask]
                    except:  # Not sure what is happening here
                        logging.info('Problems setting mask on missing_indices!')
                else:
                    missing_indices = []  # temporary workaround
                if (type(missing_indices) == np.int64) or (
                        type(missing_indices) == np.int32):
                    missing_indices = []
                self.timer_end('main loop:readers:' +
                               reader_name.replace(':', '<colon>'))
                if len(missing_indices) == 0:
                    logging.debug('Obtained data for all elements.')
                    #---------------------------------------------------------------------------------------
                    # if the current variables_group is ['x_sea_water_velocity','y_sea_water_velocity'], and 
                    # that are more than one reader available for that group , we add the velocity components
                    # rather than exiting the loop now
                    if 'x_sea_water_velocity' in variable_group and add_water_velocity :
                        pass
                    else : 
                    # standard case :    
                    # variable_group other than ['x_sea_water_velocity','y_sea_water_velocity'] should not need to be combined
                    # exit the loop now
                        break
                    #---------------------------------------------------------------------------------------
                else:
                    logging.debug('Data missing for %i elements.' %
                                  (len(missing_indices)))

        logging.debug('---------------------------------------')
        logging.debug('Finished processing all variable groups')

        self.timer_start('main loop:readers:postprocessing')
        for var in self.fallback_values:
            if (var not in variables) and (profiles is None or var not in profiles):
                continue
            mask = env[var].mask
            if sum(mask==True) > 0:
                logging.debug('    Using fallback value %s for %s for %s elements' %
                              (self.fallback_values[var], var, sum(mask==True)))
                env[var][mask] = self.fallback_values[var]
            # Profiles
            if profiles is not None and var in profiles:
                if 'env_profiles' not in locals():
                    logging.debug('Creating empty dictionary for profiles not '
                                  'profided by any reader: ' + str(self.required_profiles))
                    env_profiles = {}
                    env_profiles['z'] = \
                        np.array(self.required_profiles_z_range)[::-1]
                if var not in env_profiles:
                    logging.debug('      Using fallback value %s for %s for all profiles' %
                                  (self.fallback_values[var], var))
                    env_profiles[var] = self.fallback_values[var]*\
                        np.ma.ones((len(env_profiles['z']), self.num_elements_active()))
                else:
                    mask = env_profiles[var].mask
                    num_masked_values_per_element = sum(mask==True)
                    num_missing_profiles = sum(num_masked_values_per_element == len(env_profiles['z']))
                    env_profiles[var][mask] = self.fallback_values[var]
                    logging.debug('      Using fallback value %s for %s for %s profiles' %
                                  (self.fallback_values[var], var, num_missing_profiles,))
                    num_missing_individual = sum(num_masked_values_per_element > 0) - num_missing_profiles
                    if num_missing_individual > 0:
                        logging.debug('        ...plus %s individual points in other profiles' %
                                      num_missing_individual)

        #######################################################
        # Some extra checks of units and realistic magnitude
        #######################################################
        if 'sea_water_temperature' in variables:
            t_kelvin = np.where(env['sea_water_temperature']>100)[0]
            if len(t_kelvin) > 0:
                logging.warning('Converting temperatures from Kelvin to Celcius')
                env['sea_water_temperature'][t_kelvin] = env['sea_water_temperature'][t_kelvin] - 273.15
                if 'env_profiles' in locals() and 'sea_water_temperature' in env_profiles.keys():
                  env_profiles['sea_water_temperature'][:,t_kelvin] = \
                    env_profiles['sea_water_temperature'][:,t_kelvin] - 273.15

        #######################################################
        # Parameterisation of unavailable variables
        #######################################################
        if self.get_config('drift:use_tabularised_stokes_drift') is True:
            if 'x_wind' not in variables:
                logging.debug('No wind available to calculate Stokes drift')
            else:
                if 'sea_surface_wave_stokes_drift_x_velocity' not in variables or (
                    env['sea_surface_wave_stokes_drift_x_velocity'].max() == 0 and 
                    env['sea_surface_wave_stokes_drift_y_velocity'].max() == 0):
                        logging.info('Calculating parameterised stokes drift')
                        for i in range(len(env['x_wind'])):
                            env['sea_surface_wave_stokes_drift_x_velocity'][i], \
                            env['sea_surface_wave_stokes_drift_y_velocity'][i] = \
                                self.wave_stokes_drift_parameterised((env['x_wind'][i], env['y_wind'][i]),
                                self.get_config('drift:tabularised_stokes_drift_fetch'))

                if (env['sea_surface_wave_significant_height'].max() == 0):
                        logging.info('Calculating parameterised significant wave height')
                        for i in range(len(env['x_wind'])):
                            env['sea_surface_wave_significant_height'][i] = \
                                self.wave_significant_height_parameterised((env['x_wind'][i], env['y_wind'][i]),
                                self.get_config('drift:tabularised_stokes_drift_fetch'))
       
        #############################
        # Add uncertainty/diffusion
        #############################
        # Current
        if 'x_sea_water_velocity' in variables and \
                'y_sea_water_velocity' in variables:
            std = self.get_config('drift:current_uncertainty')
            if std > 0:
                logging.debug('Adding uncertainty for current: %s m/s' % std)
                env['x_sea_water_velocity'] += np.random.normal(
                    0, std, self.num_elements_active())
                env['y_sea_water_velocity'] += np.random.normal(
                    0, std, self.num_elements_active())
            std = self.get_config('drift:current_uncertainty_uniform')
            if std > 0:
                logging.debug('Adding uncertainty for current: %s m/s' % std)
                env['x_sea_water_velocity'] += np.random.uniform(
                    -std, std, self.num_elements_active())
                env['y_sea_water_velocity'] += np.random.uniform(
                    -std, std, self.num_elements_active())
        # Wind
        if 'x_wind' in variables and 'y_wind' in variables:
            std = self.get_config('drift:wind_uncertainty')
            if std > 0:
                logging.debug('Adding uncertainty for wind: %s m/s' % std)
                env['x_wind'] += np.random.normal(
                    0, std, self.num_elements_active())
                env['y_wind'] += np.random.normal(
                    0, std, self.num_elements_active())

        #####################
        # Diagnostic output
        #####################
        if len(env) > 0:
            logging.debug('------------ SUMMARY -------------')
            for var in variables:
                logging.debug('    %s: %g (min) %g (max)' %
                              (var, env[var].min(), env[var].max()))
            logging.debug('---------------------------------')
            logging.debug('\t\t%s active elements' % self.num_elements_active())
            if self.num_elements_active() > 0:
                lonmin = self.elements.lon.min()
                lonmax = self.elements.lon.max()
                latmin = self.elements.lat.min()
                latmax = self.elements.lat.max()
                zmin = self.elements.z.min()
                zmax = self.elements.z.max()
                if latmin == latmax:
                    logging.debug('\t\tlatitude =  %s' % (latmin))
                else:
                    logging.debug('\t\t%s <- latitude  -> %s' % (latmin, latmax))
                if lonmin == lonmax:
                    logging.debug('\t\tlongitude = %s' % (lonmin))
                else:
                    logging.debug('\t\t%s <- longitude -> %s' % (lonmin, lonmax))
                if zmin == zmax:
                    logging.debug('\t\tz = %s' % (zmin))
                else:
                    logging.debug('\t\t%s   <- z ->   %s' % (zmin, zmax))
                logging.debug('---------------------------------')

        # Prepare array indiciating which elements contain any invalid values
        missing = np.ma.masked_invalid(env[variables[0]]).mask
        for var in variables[1:]:
            missing = np.ma.mask_or(missing,
                                    np.ma.masked_invalid(env[var]).mask,
                                    shrink=False)

        # Convert dictionary to recarray and return
        if 'env_profiles' not in locals():
            env_profiles = None

        # Convert masked arrays to regular arrays for increased performance
        env = np.array(env)
        if env_profiles is not None:
            for var in env_profiles:
                env_profiles[var] = np.array(env_profiles[var])

        self.timer_end('main loop:readers:postprocessing')
        self.timer_end('main loop:readers')

        return env.view(np.recarray), env_profiles, missing



    def update_terminal_velocity(self, *args, **kwargs):
        '''
        Terminal velocity due to buoyancy or sedimentation rate,
        to be used in turbulent mixing module.
        stored as an array in self.elements.terminal_velocity
        '''
        # do nothing 
        # self.elements.terminal_velocity = 0.
        pass

    def horizontal_diffusion(self, *args, **kwargs):
        '''
        Horizontal diffusion based on random walk technique
        using diffusion coefficients K_xy (in m2/s - constant or interpolaote from field)
        and uniformly distributed random numbers R(-1,1) with a zero mean.

        
        Constant coefficients 'ocean_horizontal_diffusivity' can be set using:
        o.fallback_values['ocean_horizontal_diffusivity'] = 0.1 

        Time-Space varying coefficients  'ocean_horizontal_diffusivity' 
        can also be interpolated from an environment object (as in vertical_mixing() )
        
        * this should not be used in combination with drift:current_uncertainty ~=0 
        * as this model the same process, only with a different approach 

        The random component added to the particle position to reproduce turbulent horizontal 
        is computed using a approach similar to PyGnome, CMS :

        -https://github.com/beatrixparis/connectivity-modeling-system/blob/master/cms-master/src/mod_turb.f90
        -Garcia-Martinez and Tovar, 1999 - Computer Modeling of Oil Spill Trajectories With a High Accuracy Method
        -Lonin, S.A., 1999. Lagrangian model for oil spill diffusion at sea. Spill Science and Technology Bulletin, 5(5): 331-336
        -https://github.com/NOAA-ORR-ERD/PyGnome/blob/master/lib_gnome/Random_c.cpp#L50 

        stored as an array in self.elements.horizontal_diffusion
 
        '''
        if not self.environment.ocean_horizontal_diffusivity.any():
            logging.debug('No horizontal diffusion applied - ocean_horizontal_diffusivity = 0.0')
            pass

        # check if some diffusion is not already accounted for using drift:current_uncertainty
        if self.get_config('drift:current_uncertainty') != 0:
            logging.debug('Warning = some horizontal diffusion already accounted for using drift:current_uncertainty')

        diff_fac = 6 # 6 is used in PyGnome as well, but can be = 2 in other formulations, such as in the CMS model - hard coded for now
        K_xy = self.environment.ocean_horizontal_diffusivity # horizontal diffusion coefficient in [m2/s]
        # max diffusion distance in meters is : max_diff_distance = np.sqrt(diff_fac * K_xy * self.time_step.seconds)
        # here we need velocity to fit with the update_position() subroutine
        # max_diff_velocity = max_diff_distance / self.time_step.seconds = np.sqrt(diff_fac * K_xy / self.time_step.seconds) 
        max_diff_velocity = np.sqrt(diff_fac * K_xy / self.time_step.seconds) # max diffusion velocity in [m/s]
        # add randomness - random number in [-1,1] using uniform distribution i.e. same probabilty for each value (note some other formulations may use "normal" distribution)
        diff_velocity = np.random.uniform(-1,1,size = len(max_diff_velocity)) * max_diff_velocity # in [m/s]
        # random directions
        theta_rand = np.random.uniform( 0, 2*np.pi, size = len(max_diff_velocity) )
        # split diff_velocity into (u,v) velocities using random directions
        x_vel = diff_velocity * np.cos(theta_rand)
        y_vel = diff_velocity * np.sin(theta_rand)

        logging.debug('Applying horizontal diffusion to particle positions')
        logging.debug('\t\t%s   <- horizontal diffusion distance [m] ->   %s' % (np.min(diff_velocity*self.time_step.seconds), np.max(diff_velocity*self.time_step.seconds)))

        # update positions with the diffusion velocities     
        self.update_positions(x_vel, y_vel)

    def update(self):
        """Update positions and properties of elements."""

        self.elements.age_seconds += self.time_step.total_seconds()

        # Simply move particles with ambient current
        self.advect_ocean_current() # from physics_methods.py
        # Horizontal diffusion
        self.horizontal_diffusion()

        # Advect particles due to wind drag
        # (according to specified wind_drift_factor)
        self.advect_wind()

        # Stokes drift
        self.stokes_drift()

        # Turbulent Mixing
        self.update_terminal_velocity()
        self.vertical_mixing() # using subroutine from opendrift3D.py

        # Vertical advection
        self.vertical_advection()

        # Sediment resuspension checks , if switched on
        # self.sediment_resuspension() - ToDo!
        # 
        # 1-find particles on the bottom
        # 2-compute bed shear stresses
        # 3-compare to critical_shear_stress
        # 4-resuspend or stay on seabed depending on 3)
        #   > probably need to use a cut-off age after which particles are de-activated anyway
        #   to prevent excessive build-up of "active" particle in the simulations

        # Deactivate elements that exceed a certain age
        if self.get_config('drift:max_age_seconds') is not None:
            self.deactivate_elements(self.elements.age_seconds >=
                                     self.get_config('drift:max_age_seconds'),
                                     reason='retired')
        # When no resuspension is required, deactivate that reached the seabed
        if self.get_config('processes:resuspension') is False:
            self.deactivate_elements(self.elements.z ==
                                     -1.*self.environment.sea_floor_depth_below_sea_level,
                                     reason='settled')

        # Note the interaction with shoreline in taken care of by interact_with_coastline in basemodel.py
        # when run() is called
