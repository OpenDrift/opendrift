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
# Copyright 2020, Gaute Hope, MET Norway

import logging
logging.captureWarnings(True)
logger = logging.getLogger('opendrift')
logging.getLogger('botocore').setLevel(logging.INFO)
logging.getLogger('urllib3').setLevel(logging.INFO)

import sys
import os
import types
from typing import Union, List
import traceback
import inspect
import psutil

from opendrift.models.basemodel.environment import Environment
from opendrift.readers import reader_global_landmask

from datetime import datetime, timedelta
from abc import ABCMeta, abstractmethod, abstractproperty

import geojson
import xarray as xr
import pandas as pd
import numpy as np
import scipy
import pyproj
import matplotlib

matplotlib.pyplot.set_loglevel(level = 'warning')
matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['legend.scatterpoints'] = 1
matplotlib.rcParams['figure.autolayout'] = True
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Polygon
from matplotlib.path import Path
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from enum import Enum
import functools

import opendrift
from opendrift.timer import Timeable
from opendrift.errors import WrongMode
from opendrift.models.physics_methods import PhysicsMethods
from opendrift.config import Configurable, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED, CONFIG_LEVEL_ESSENTIAL

import roaring_landmask
from roaring_landmask import RoaringLandmask

Mode = Enum('Mode', ['Config', 'Ready', 'Run', 'Result'])
rl = roaring_landmask.RoaringLandmask.new()

def require_mode(mode: Union[Mode, List[Mode]], post_next_mode=False, error=None):
    if not isinstance(mode, list):
        mode = [mode]

    def _decorator(func):

        @functools.wraps(func)
        def inner(self, *args, **kwargs):

            def next_mode():
                # Change the mode
                prev = self.mode

                if self.mode is Mode.Config:
                    self.env.finalize()
                    self.mode = Mode.Ready

                elif self.mode is Mode.Ready:
                    self.mode = Mode.Run

                elif self.mode is Mode.Run:
                    self.mode = Mode.Result

                elif self.mode is Mode.Result:
                    pass

                else:
                    raise Exception("Unknown mode")

                logger.debug(f"Changed mode from {prev} to {self.mode}")

            if self.mode not in mode:
                # Check if we can advance to the required mode
                if mode[0] is Mode.Ready and self.mode is Mode.Config:
                    next_mode()

                elif mode[0] is Mode.Run and self.mode is Mode.Ready:
                    next_mode()

                elif mode[0] is Mode.Result and self.mode is Mode.Run:
                    next_mode()

                else:
                    raise WrongMode(mode, self.mode, error)

            r = func(self, *args, **kwargs)

            if post_next_mode:
                next_mode()

            return r

        return inner

    return _decorator


class OpenDriftSimulation(PhysicsMethods, Timeable, Configurable):
    """Generic trajectory model class, to be extended (subclassed).

    This as an Abstract Base Class, meaning that only subclasses can
    be initiated and used.
    Any specific subclass ('model') must contain its own (or shared)
    specific type of particles (ElementType), whose properties are
    updated at each time_step using method update() on basis of model
    physics/chemistry/biology and 'required_variables' (environment)
    which are provided by one or more Reader objects.

    Attributes:
        ElementType: the type (class) of particles to be used by this model

        elements: object of the class ElementType, storing the specific
        particle properties (ndarrays and scalars) of all active particles
        as named attributes. Elements are added by seeding-functions
        (presently only one implemented: seed_elements).

        elements_deactivated: ElementType object containing particles which
            have been deactivated (and removed from 'elements')

        elements_scheduled: ElementType object containing particles which
            have been scheduled, but not yet activated
        required_variables: list of strings of CF standard_names which is
            needed by this model (update function) to update properties of
            particles ('elements') at each time_step. This core class has
            no required_elements, this is implemented by subclasses/modules.
        environment: recarray storing environment variables (wind, waves,
            current etc) as named attributes. Attribute names follow
            standard_name from CF-convention, allowing any OpenDriftSimulation
            module/subclass using environment data from any readers which
            can provide the requested variables. Used in method 'update'
            to update properties of elements every time_step.
        time_step: timedelta object, time interval at which element properties
            are updated (including advection).
        time_step_output: timedelta object, time interval at which element
            properties are stored in memory and eventually written to file
        readers: Dictionary where values are Reader objects, and names are
            unique reference keywords used to access a given reader (typically
            filename or URL)
        priority_list: OrderedDict where names are variable names,
            and values are lists of names (kewywords) of the reader, in the
            order of priority (user defined) of which the readers shall be
            called to retrieve environmental data.
    """

    __metaclass__ = ABCMeta

    mode = Mode.Config

    status_categories = ['active']  # Particles are active by default

    # Default plotting colors of trajectory endpoints
    status_colors_default = {
        'initial': 'green',
        'active': 'blue',
        'missing_data': 'gray'
    }

    plot_comparison_colors = [
        'k', 'r', 'g', 'b', 'm', 'c', 'y', 'crimson', 'indigo', 'lightcoral',
        'grey', 'sandybrown', 'palegreen', 'gold', 'yellowgreen', 'lime',
        'steelblue', 'navy', 'darkviolet'
    ]
    plot_comparison_colors = plot_comparison_colors + plot_comparison_colors

    proj_latlon = pyproj.Proj('+proj=latlong')
    """
    The environment holds all the readers and the forcing data for the simulation.
    """
    env: Environment

    @classmethod
    def SRS(cls):
        return cls.proj_latlon

    def __init__(self,
                 seed=0,
                 iomodule='netcdf',
                 loglevel=logging.DEBUG,
                 logfile=None):
        """Initialise OpenDriftSimulation

        Args:
            seed: integer or None. A given integer will yield identical
                random numbers drawn each simulation. Random numbers are
                e.g. used to distribute particles spatially when seeding,
                and may be used by modules (subclasses) for e.g. diffusion.
                Specifying a fixed value (default: 0) is useful for sensitivity
                tests. With seed = None, different random numbers will be drawn
                for subsequent runs, even with identical configuration/input.
            iomodule: name of module used to export data
                default: netcdf, see :py:mod:`opendrift.io` for more alternatives.
                `iomodule` is module/filename without preceeding `io_`
            loglevel: set to 0 (default) to retrieve all debug information.
                Provide a higher value (e.g. 20) to receive less output.
            logtime: if True, a time stamp is given for each logging line.
                     logtime can also be given as a python time specifier
                     (e.g. '%H:%M:%S')
            logfile:
                None (default) to send output to console.
                A string to send output to logfile.
                Or to get output to both terminal and file:
                    [<filename.log>, logging.StreamHandler(sys.stdout)]
        """

        super().__init__()

        self.profiles_depth = None

        self.show_continuous_performance = False

        self.origin_marker = None  # Dictionary to store named seeding locations

        # List to store GeoJSON dicts of seeding commands
        self.seed_geojson = []

        self.env = Environment(self.required_variables, self._config)

        # Make copies of dictionaries so that they are private to each instance
        self.status_categories = ['active']  # Particles are active by default
        self.status_colors_default = self.status_colors_default.copy()

        if hasattr(self, 'status_colors'):
            # Append model specific colors to (and override) default colors
            self.status_colors_default.update(self.status_colors)
            self.status_colors = self.status_colors_default
        else:
            self.status_colors = self.status_colors_default

        # Using a fixed seed will generate the same random numbers
        # each run, useful for sensitivity tests
        # Use seed = None to get different random numbers each time
        np.random.seed(seed)

        self.steps_calculation = 0  # Increase for each simulation step
        self.elements_deactivated = self.ElementType()  # Empty array
        self.elements = self.ElementType()  # Empty array

        # Set up logging
        logformat = '%(asctime)s %(levelname)-7s %(name)s:%(lineno)d: %(message)s'
        datefmt = '%H:%M:%S'

        if loglevel < 10:  # 0 is NOTSET, giving no output
            loglevel = 10
        logger.setLevel(loglevel)
        logger.handlers.clear()

        if logfile is not None:
            if not isinstance(logfile, list):
                logfile = [logfile]
            for lh in logfile:
                if not isinstance(lh, logging.StreamHandler):
                    lh = logging.FileHandler(lh, mode='w')
                lh.setLevel(loglevel)
                lh.setFormatter(logging.Formatter(logformat, datefmt))
                logger.addHandler(lh)
        else:  # logging only to console, with colors
            import coloredlogs
            fields = coloredlogs.DEFAULT_FIELD_STYLES
            fields['levelname']['color'] = 'magenta'

            # coloredlogs does not create duplicate handlers
            coloredlogs.install(level=loglevel,
                                fmt=logformat,
                                datefmt=datefmt,
                                field_styles=fields)

        # Prepare outfile
        try:
            io_module = __import__(
                'opendrift.export.io_' + iomodule,
                fromlist=['init', 'write_buffer', 'close', 'import_file'])
        except ImportError:
            logger.info('Could not import iomodule ' + iomodule)
            raise
        self.io_init = types.MethodType(io_module.init, self)
        self.io_write_buffer = types.MethodType(io_module.write_buffer, self)
        self.io_close = types.MethodType(io_module.close, self)
        self.io_import_file = types.MethodType(io_module.import_file, self)

        # Set configuration options
        self._add_config({
            # type, default, min, max, enum, important, value, units, description
            'general:simulation_name': {'type': 'str', 'min_length': 0, 'max_length': 64,
                                        'default': '', 'level': CONFIG_LEVEL_BASIC,
                                        'description': 'Name of simulation'},
            'general:coastline_action': {
                'type': 'enum',
                'enum': ['none', 'stranding', 'previous'],
                'default': 'stranding',
                'level': CONFIG_LEVEL_BASIC,
                'description': 'None means that objects may also move over land. '
                    'stranding means that objects are deactivated if they hit land. '
                    'previous means that objects will move back to the previous location '
                    'if they hit land'
            },
            'general:coastline_approximation_precision': {
                'type': 'float',
                'default': 0.001,
                'min': 0.0001,
                'max': 0.005,
                'units': 'degrees',
                'description': 'The precision of the particle position approximation to the coastline.',
                'level': CONFIG_LEVEL_BASIC
            },
            'general:time_step_minutes': {
                'type': 'float',
                'min': .01,
                'max': 1440,
                'default': 60,
                'units': 'minutes',
                'level': CONFIG_LEVEL_BASIC,
                'description': 'Calculation time step used for the simulation. The output time step may '
                'be equal or larger than this.'
            },
            'general:time_step_output_minutes': {
                'type': 'float',
                'min': 1,
                'max': 1440,
                'default': None,
                'units': 'minutes',
                'level': CONFIG_LEVEL_BASIC,
                'description': 'Output time step, i.e. the interval at which output is saved. '
                'This must be larger than the calculation time step, and be an integer multiple of this.'
            },
            'seed:ocean_only': {
                'type': 'bool',
                'default': True,
                'description': 'If True, elements seeded on land will be moved to the closest '
                'position in ocean',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'seed:number': {
                'type': 'int',
                'default': 1,
                'min': 1,
                'max': 100000000,
                'units': 1,
                'description': 'The number of elements for the simulation.',
                'level': CONFIG_LEVEL_BASIC
            },
            'drift:max_age_seconds': {
                'type': 'float',
                'default': None,
                'min': 0,
                'max': np.inf,
                'units': 'seconds',
                'description':
                'Elements will be deactivated when this age is reached',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:advection_scheme': {
                'type': 'enum',
                'enum': ['euler', 'runge-kutta', 'runge-kutta4'],
                'default': 'euler',
                'level': CONFIG_LEVEL_ADVANCED,
                'description': 'Numerical advection scheme for ocean current advection'
            },
            'drift:horizontal_diffusivity': {
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 100000,
                'units': 'm2/s',
                'description': 'Add horizontal diffusivity (random walk)',
                'level': CONFIG_LEVEL_BASIC
            },
            'drift:profiles_depth': {'type': 'float', 'default': 50, 'min': 0, 'max': None,
                'level': CONFIG_LEVEL_ADVANCED, 'units': 'meters', 'description':
                'Environment profiles will be retrieved from surface and down to this depth'},
            'drift:wind_uncertainty': {
                'type': 'float',
                'default': 0,
                'min': 0,
                'max': 5,
                'units': 'm/s',
                'description':
                'Add gaussian perturbation with this standard deviation to wind components at each time step.',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:relative_wind': {
                'type': 'bool',
                'default': False,
                'description':
                'If True, wind drift is calculated for absolute wind (wind vector minus ocean surface current vector).',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:deactivate_north_of': {
                'type': 'float',
                'default': None,
                'min': -90,
                'max': 90,
                'units': 'degrees',
                'description':
                'Elements are deactivated if the move further north than this limit',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:deactivate_south_of': {
                'type': 'float',
                'default': None,
                'min': -90,
                'max': 90,
                'units': 'degrees',
                'description':
                'Elements are deactivated if the move further south than this limit',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:deactivate_east_of': {
                'type': 'float',
                'default': None,
                'min': -360,
                'max': 360,
                'units': 'degrees',
                'description':
                'Elements are deactivated if the move further east than this limit',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'drift:deactivate_west_of': {
                'type': 'float',
                'default': None,
                'min': -360,
                'max': 360,
                'units': 'degrees',
                'description':
                'Elements are deactivated if the move further west than this limit',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'readers:max_number_of_fails': {
                'type': 'int',
                'default': 1,
                'min': 0,
                'max': 1e6,
                'units': 'number',
                'description':
                'Readers are discarded if they fail (e.g. corrupted data, og hanging servers) move than this number of times',
                'level': CONFIG_LEVEL_ADVANCED
            },
        })

        # Add default element properties to config
        c = {}
        for p in self.ElementType.variables:
            v = self.ElementType.variables[p]
            if 'seed' in v and v['seed'] is False:
                continue  # Properties which may not be provided by user
            c['seed:%s' % p] = {
                'type': v['type'] if 'type' in v else 'float',
                'min': v['min'] if 'min' in v else None,
                'max': v['max'] if 'max' in v else None,
                'units': v['units'] if 'units' in v else None,
                'default': v['default'] if 'default' in v else None,
                'description': v['description'] if 'description' in v \
                        else 'Seeding value of %s' % p, 'level': v['level'] if 'level' in v \
                        else CONFIG_LEVEL_ADVANCED
            }
        self._add_config(c)

        self.result = None  # Xarray Dataset to store trajectories and properties
        self.index_of_first = None
        self.index_of_last = None

        # Find variables which require profiles
        self.required_profiles = [
            var for var in self.required_variables
            if 'profiles' in self.required_variables[var]
            and self.required_variables[var]['profiles'] is True
        ]

        # Find variables which are desired, but not required
        self.desired_variables = [
            var for var in self.required_variables
            if 'important' in self.required_variables[var]
            and self.required_variables[var]['important'] is False
        ]

        self.timer_start('total time')
        self.timer_start('configuration')

        self.add_metadata('opendrift_version', opendrift.__version__)
        logger.info('OpenDriftSimulation initialised (version %s)' %
                    opendrift.version.version_or_git())

    def clone(self):
        c = self.__class__()

        c._config.clear()
        for k, v in self._config.items():
            c._config[k] = v

        c.add_reader([r for _, r in self.env.readers.items()])

        return c

    @require_mode(mode=Mode.Config,
                  error='Cannot set config after elements have been seeded')
    @functools.wraps(Configurable.set_config)
    def set_config(self, *args, **kwargs):
        return Configurable.set_config(self, *args, **kwargs)

    @require_mode(mode=[Mode.Config, Mode.Ready])
    @functools.wraps(Configurable.set_config)
    def __set_seed_config__(self, key: str, value):
        """
        This method allows setting config values that are passed as seed arguments. The environment is already prepared before this, so make sure that nothing is changed that requires the environment to be re-initialized.
        """
        if not key.startswith('seed'):
            raise ValueError("This method is only allowed for setting seed arguments.")

        # check that the oil_type is only set once
        if key == 'seed:oil_type' and self.num_elements_scheduled() > 0:
            if value != self.get_config('seed:oil_type'):
                raise WrongMode(Mode.Config, self.mode, msg=f"Cannot change the oil type after elements have been seeded: {self.get_config('seed:oil_type')} -> {value}")

        return Configurable.set_config(self, key, value)

    def add_metadata(self, key, value):
        """Add item to metadata dictionary, for export as netCDF global attributes"""
        if not hasattr(self, 'metadata_dict'):
            from collections import OrderedDict
            self.metadata_dict = OrderedDict()
        self.metadata_dict[key] = value

    @require_mode(mode=[Mode.Config, Mode.Result])
    def add_reader(self, readers, variables=None, first=False):
        self.env.add_reader(readers, variables, first)

    @require_mode(mode=Mode.Config)
    def add_readers_from_list(self, *args, **kwargs):
        '''Make readers from a list of URLs or paths to netCDF datasets'''
        self.env.add_readers_from_list(*args, **kwargs)

    @require_mode(mode=Mode.Config)
    def add_readers_from_file(self, *args, **kwargs):
        '''Make readers from a file containing list of URLs or paths to netCDF datasets'''
        self.env.add_readers_from_file(*args, **kwargs)

    # To be overloaded by sublasses, but this parent method must be called
    def prepare_run(self):
        # Copy profile_depth from config
        self.profiles_depth = self.get_config('drift:profiles_depth')

    # To be overloaded by sublasses, but this parent method must be called
    def post_run(self):
        pass

    def store_present_positions(self, IDs=None, lons=None, lats=None):
        """Store present element positions, in case they shall be moved back"""
        if self.get_config('general:coastline_action') in ['previous', 'stranding'] or (
                'general:seafloor_action' in self._config
                and self.get_config('general:seafloor_action') == 'previous'):
            if not hasattr(self, 'previous_lon'):
                self.previous_lon = np.ma.masked_all(self.num_elements_total())
                self.previous_lat = np.ma.masked_all(self.num_elements_total())
            if IDs is None:
                IDs = self.elements.ID
                lons = self.elements.lon
                lats = self.elements.lat
                self.newly_seeded_IDs = None
            else:
                # to check if seeded on land
                if len(IDs) > 0:
                    self.newly_seeded_IDs = np.copy(IDs)
                else:
                    self.newly_seeded_IDs = None
            self.previous_lon[IDs] = np.copy(lons)
            self.previous_lat[IDs] = np.copy(lats)

    def store_previous_variables(self):
        """Store some environment variables, for access at next time step"""

        if not hasattr(self, 'store_previous'):
            return
        if not hasattr(self, 'variables_previous'):
            # Create ndarray to store previous variables
            dtype = [(var, np.float32) for var in self.store_previous]
            self.variables_previous = np.array(np.full(
                self.num_elements_total(), np.nan),
                                               dtype=dtype)

        # Copying variables_previous to environment_previous
        self.environment_previous = self.variables_previous[self.elements.ID]

        # Use new values for new elements which have no previous value
        for var in self.store_previous:
            undefined = np.isnan(self.environment_previous[var])
            self.environment_previous[var][undefined] = getattr(
                self.environment, var)[undefined]

        self.environment_previous = self.environment_previous.view(np.recarray)

        for var in self.store_previous:
            self.variables_previous[var][self.elements.ID] = getattr(
                self.environment, var)

    def interact_with_coastline(self, final=False):
        """Coastline interaction according to configuration setting"""
        if self.num_elements_active() == 0:
            return
        i = self.get_config('general:coastline_action')
        coastline_approximation_precision = self.get_config('general:coastline_approximation_precision')
        if not hasattr(self, 'environment') or not hasattr(
                self.environment, 'land_binary_mask'):
            return
        if i == 'none':  # Do nothing
            return
        if final is True:  # Get land_binary_mask for final location
            en, en_prof, missing = \
                self.env.get_environment(['land_binary_mask'],
                                     self.time,
                                     self.elements.lon,
                                     self.elements.lat,
                                     self.elements.z)
            self.environment.land_binary_mask = en.land_binary_mask

        if i == 'stranding':  # Deactivate elements on land, but not in air
            on_land = np.where(self.environment.land_binary_mask == 1)[0]
            if len(on_land) == 0:
                logger.debug('No elements hit coastline.')
                return

            logger.debug('%s elements hit land, moving them to the coastline.' % len(on_land))

            self.deactivate_elements(
                (self.environment.land_binary_mask == 1) & (self.elements.z <= 0),
                reason='stranded'
            )

            if not coastline_approximation_precision:
                return

            for on_land_id, on_land_prev_id in zip(on_land, self.elements.ID[on_land]):
                lon = self.elements.lon[on_land_id]
                lat = self.elements.lat[on_land_id]
                prev_lon = self.previous_lon[on_land_prev_id]
                prev_lat = self.previous_lat[on_land_prev_id]

                step_degrees = float(coastline_approximation_precision)

                x_degree_diff = np.abs(prev_lon - lon)
                x_samples = np.floor(x_degree_diff / step_degrees).astype(np.int64) if x_degree_diff > step_degrees else 1
                x = np.linspace(prev_lon, lon, x_samples)

                y_degree_diff = np.abs(prev_lat - lat)
                y_samples = np.floor(y_degree_diff/ step_degrees).astype(np.int64) if y_degree_diff > step_degrees else 1
                y = np.linspace(prev_lat, lat, y_samples)

                xx, yy = np.meshgrid(x,y)
                xx, yy = xx.ravel(), yy.ravel()

                rl_mask = rl.contains_many(xx.ravel(), yy.ravel())
                if np.any(rl_mask):
                    index = np.argmax(rl_mask)
                    new_lon = xx[index]
                    new_lat = yy[index]

                    self.elements.lon[on_land_id] = new_lon
                    self.elements.lat[on_land_id] = new_lat

            self.environment.land_binary_mask[on_land] = 0

        elif i == 'previous':  # Go back to previous position (in water)
            if self.newly_seeded_IDs is not None:
                self.deactivate_elements(
                    (self.environment.land_binary_mask == 1) &
                    (self.elements.age_seconds == 0),
                    reason='seeded_on_land')
            on_land = np.where(self.environment.land_binary_mask == 1)[0]
            if len(on_land) == 0:
                logger.debug('No elements hit coastline.')
            else:
                logger.debug('%s elements hit coastline, '
                             'moving back to water' % len(on_land))
                on_land_ID = self.elements.ID[on_land]
                self.elements.lon[on_land] = \
                    np.copy(self.previous_lon[on_land_ID])
                self.elements.lat[on_land] = \
                    np.copy(self.previous_lat[on_land_ID])
                self.environment.land_binary_mask[on_land] = 0

    def interact_with_seafloor(self):
        """Seafloor interaction according to configuration setting"""
        if self.num_elements_active() == 0:
            return
        if 'sea_floor_depth_below_sea_level' not in self.env.priority_list:
            return

        # The shape of these is different than the original arrays
        # because it is for active drifters
        sea_floor_depth = self.sea_floor_depth()
        sea_surface_height = self.sea_surface_height()

        # Check if any elements are below sea floor
        # But remember that the water column is the sea floor depth + sea surface height
        ibelow = self.elements.z < -(sea_floor_depth + sea_surface_height)
        below = np.where(ibelow)[0]

        if len(below) == 0:
            logger.debug('No elements hit seafloor.')
            return

        i = self.get_config('general:seafloor_action')
        if i == 'lift_to_seafloor':
            logger.debug('Lifting %s elements to seafloor.' % len(below))
            self.elements.z[below] = -(sea_floor_depth + sea_surface_height)[below]
        elif i == 'deactivate':
            self.deactivate_elements(ibelow, reason='seafloor')
            self.elements.z[below] = -(sea_floor_depth + sea_surface_height)[below]
        elif i == 'previous':  # Go back to previous position (in water)
            logger.warning('%s elements hit seafloor, '
                           'moving back ' % len(below))
            below_ID = self.elements.ID[below]
            self.elements.lon[below] = \
                np.copy(self.previous_lon[below_ID])
            self.elements.lat[below] = \
                np.copy(self.previous_lat[below_ID])

    @abstractmethod
    def update(self):
        """Any trajectory model implementation must define an update method.
        This method must/can use environment data (self.environment) to
        update properties (including position) of its particles (self.elements)
        """

    @abstractproperty
    def ElementType(self):
        """Any trajectory model implementation must define an ElementType."""

    @abstractproperty
    def required_variables(self):
        """Any trajectory model implementation must list needed variables."""

    def test_data_folder(self):
        import opendrift
        return os.path.abspath(
            os.path.join(os.path.dirname(opendrift.__file__), '..', 'tests',
                         'test_data')) + os.path.sep

    def performance(self):
        '''Report the time spent on various tasks'''

        outStr = '--------------------\n'
        outStr += 'Reader performance:\n'
        for r in self.env.readers:
            reader = self.env.readers[r]
            if reader.is_lazy:
                continue
            outStr += '--------------------\n'
            outStr += r + '\n'
            outStr += reader.performance()

        outStr += '--------------------\n'
        outStr += 'Performance:\n'
        for category, time in self.timing.items():
            timestr = str(time)[0:str(time).find('.') + 2]
            for i, c in enumerate(timestr):
                if c in '123456789.':
                    timestr = timestr[i:]  # Strip leading 0 and :
                    if c == '.':
                        timestr = '0' + timestr
                    break
            parts = category.split(':')
            indent = '  ' * (len(parts) - 1)
            category = parts[-1]
            category = category.replace('<colon>', ':')
            outStr += '%s%7s %s\n' % (indent, timestr, category)

        outStr += '--------------------\n'
        return outStr

    def num_elements_active(self):
        """The number of active elements."""
        if hasattr(self, 'elements'):
            return len(self.elements)
        else:
            return 0

    def num_elements_deactivated(self):
        """The number of deactivated elements."""
        if hasattr(self, 'elements_deactivated'):
            return len(self.elements_deactivated)
        else:
            return 0

    def num_elements_scheduled(self):
        if hasattr(self, 'elements_scheduled'):
            return len(self.elements_scheduled)
        else:
            return 0

    def num_elements_total(self):
        """The total number of scheduled, active and deactivated elements."""
        return self.num_elements_activated() + self.num_elements_scheduled()

    def num_elements_activated(self):
        """The total number of active and deactivated elements."""
        return self.num_elements_active() + self.num_elements_deactivated()

    @require_mode(mode=Mode.Ready)
    def schedule_elements(self, elements, time):
        """Schedule elements to be seeded during runtime.

        Also assigns a unique ID to each particle, monotonically increasing."""

        # prepare time
        if isinstance(time, np.ndarray):
            time = list(time)
        if not isinstance(time, list):
            time = [time]
        if len(time) == 1 and len(elements) > 1:
            time = time * len(elements)

        if not hasattr(self, 'elements_scheduled'):
            self.elements_scheduled = elements
            self.elements_scheduled_time = np.array(time)
            # We start simulation at time of release of first element:
            self.start_time = time[0]
            self.elements_scheduled.ID = np.arange(0, len(elements))
        else:
            elements.ID = np.arange(self.num_elements_scheduled(),
                                    self.num_elements_scheduled() +
                                    len(elements))  # Increase ID successively
            self.elements_scheduled.extend(elements)
            self.elements_scheduled_time = np.append(
                self.elements_scheduled_time, np.array(time))

        min_time = np.min(time)
        if hasattr(self, 'start_time'):
            if min_time < self.start_time:
                self.start_time = min_time
                logger.debug('Setting simulation start time to %s' %
                             str(min_time))
        else:
            self.start_time = min_time
            logger.debug('Setting simulation start time to %s' % str(min_time))

    def release_elements(self):
        """Activate elements which are scheduled within following timestep."""

        logger.debug(
            'to be seeded: %s, already seeded %s' %
            (len(self.elements_scheduled), self.num_elements_activated()))
        if len(self.elements_scheduled) == 0:
            return
        if self.time_step.days >= 0:
            indices = (self.elements_scheduled_time >= self.time) & \
                      (self.elements_scheduled_time <
                       self.time + self.time_step)
        else:
            indices = (self.elements_scheduled_time <= self.time) & \
                      (self.elements_scheduled_time >
                       self.time + self.time_step)
        self.store_present_positions(self.elements_scheduled.ID[indices],
                                     self.elements_scheduled.lon[indices],
                                     self.elements_scheduled.lat[indices])
        self.elements_scheduled.move_elements(self.elements, indices)
        self.elements_scheduled_time = self.elements_scheduled_time[~indices]
        logger.debug('Released %i new elements.' % np.sum(indices))

    def closest_ocean_points(self, lon, lat):
        """Return the closest ocean points for given lon, lat"""
        from opendrift.readers.reader_shape import Reader as ShapeReader

        if not 'land_binary_mask' in self.env.priority_list:
            logger.info('No land reader added, '
                        'making a temporary landmask reader')
            from opendrift.models.oceandrift import OceanDrift
            reader_landmask = reader_global_landmask.Reader()
            seed_state = np.random.get_state(
            )  # Do not alter current random number generator
            o = OceanDrift()
            np.random.set_state(seed_state)
            if hasattr(self, 'simulation_extent'):
                o.simulation_extent = self.simulation_extent
            o.env.add_reader(reader_landmask)
            o.env.finalize()  # This is not env of the main simulation
            land_reader = reader_landmask
        else:
            logger.info('Using existing reader for land_binary_mask')
            land_reader_name = self.env.priority_list['land_binary_mask'][0]
            land_reader = self.env.readers[land_reader_name]
            o = self

        if isinstance(land_reader, ShapeReader):
            # can do this better
            is_inside = land_reader.__on_land__(lon, lat)
            lon[is_inside], lat[is_inside], _ = land_reader.get_nearest_outside(
                lon[is_inside],
                lat[is_inside],
                buffer_distance=land_reader._land_kdtree_buffer_distance,
            )

        else:
            deltalon = 0.01  # grid
            deltalat = 0.01
            numbuffer = 10
            lonmin = lon.min() - deltalon * numbuffer
            lonmax = lon.max() + deltalon * numbuffer
            latmin = lat.min() - deltalat * numbuffer
            latmax = lat.max() + deltalat * numbuffer

            land = o.env.get_environment(['land_binary_mask'],
                                        lon=lon,
                                        lat=lat,
                                        z=0 * lon,
                                        time=land_reader.start_time)[0]['land_binary_mask']
            if land.max() == 0:
                logger.info('All points are in ocean')
                return lon, lat
            logger.info('Moving %i out of %i points from land to water' %
                        (np.sum(land != 0), len(lon)))
            landlons = lon[land != 0]
            landlats = lat[land != 0]
            longrid = np.arange(lonmin, lonmax, deltalon)
            latgrid = np.arange(latmin, latmax, deltalat)
            longrid, latgrid = np.meshgrid(longrid, latgrid)
            longrid = longrid.ravel()
            latgrid = latgrid.ravel()
            # Remove grid-points not covered by this reader
            latgrid_covered = land_reader.covers_positions(longrid, latgrid)[0]
            longrid = longrid[latgrid_covered]
            latgrid = latgrid[latgrid_covered]
            landgrid = o.env.get_environment(['land_binary_mask'],
                                            lon=longrid,
                                            lat=latgrid,
                                            z=0 * longrid,
                                            time=land_reader.start_time)[0]['land_binary_mask']
            if landgrid.size == 0:
                # Need to catch this before trying .min() on it...
                logger.warning('Land grid has zero size, cannot move elements.')
                return lon, lat                                        
            if landgrid.min() == 1 or np.isnan(landgrid.min()):
                logger.warning('No ocean pixels nearby, cannot move elements.')
                return lon, lat

            oceangridlons = longrid[landgrid == 0]
            oceangridlats = latgrid[landgrid == 0]
            tree = scipy.spatial.cKDTree(
                np.dstack([oceangridlons, oceangridlats])[0])
            landpoints = np.dstack([landlons, landlats])
            _dist, indices = tree.query(landpoints)
            indices = indices.ravel()
            lon[land != 0] = oceangridlons[indices]
            lat[land != 0] = oceangridlats[indices]

        return lon, lat

    @require_mode(mode=Mode.Ready)
    def seed_elements(self,
                      lon,
                      lat,
                      time,
                      radius=0,
                      number=None,
                      number_per_point=None,
                      radius_type='gaussian',
                      **kwargs):
        """Seed elements with given position(s), time and properties.

        Arguments:
            lon: scalar or array
                central longitude(s).
            lat: scalar or array
                central latitude(s).
            radius: scalar or array
                radius in meters around each lon-lat pair,
                within which particles will be randomly seeded.
            number: integer, total number of particles to be seeded
                If number is None, the number of elements is the
                length of lon/lat or time if these are arrays. Otherwise
                the number of elements are obtained from the config-default.
                If provided, number must be a multiple of the number of points.
            number_per_point: integer, number of particles to be seeded at each point.
                This shall not be provided along with number. 
                Only relevant if lon/lat are arrays.
            time: datenum or list
                The time at which particles are seeded/released.
                If time is a list with two elements, elements are seeded
                continously from start/first to end/last time.
                If time is a list with more than two elements, the number of elements
                is equal to len(time) and are seeded as a time series.
            radius_type: string
                If 'gaussian' (default), the radius is the standard deviation in
                x-y-directions. If 'uniform', elements are spread evenly and
                always inside a circle with the given radius.
            kwargs:
                keyword arguments containing properties/attributes and
                values corresponding to the actual particle type (ElementType).
                These are forwarded to the ElementType class. All properties
                for which there are no default value must be specified.
        """

        if 'cone' in kwargs:
            raise ValueError(
                'Keyword *cone* for seed_elements is deprecated, use seed_cone() instead.'
            )

        if self.origin_marker is None:
            self.origin_marker = {}
        if 'origin_marker' in kwargs:
            origin_marker = kwargs['origin_marker']
        else:
            origin_marker = len(self.origin_marker)
        if 'origin_marker_name' in kwargs:
            origin_marker_name = kwargs['origin_marker_name']
            del kwargs['origin_marker_name']
        else:
            origin_marker_name = 'Seed %d' % len(self.origin_marker)
        if not 'origin_marker' in kwargs:
            kwargs['origin_marker'] = origin_marker
        if '_' in origin_marker_name:
            raise ValueError(
                'Underscore (_) not allowed in origin_marker_name')
        self.origin_marker[str(origin_marker)] = origin_marker_name.replace(
            ' ', '_')

        lon = np.atleast_1d(lon).ravel()
        lat = np.atleast_1d(lat).ravel()
        radius = np.atleast_1d(radius).ravel()
        time = np.atleast_1d(time)

        if lat.max() > 90 or lat.min() < -90:
            raise ValueError('Latitude must be between -90 and 90 degrees')

        if len(lon) != len(lat):
            raise ValueError('Lon and lat must have same lengths')

        if len(lon) > 1:
            if number_per_point is not None:
                if number is not None:
                    raise ValueError('Both number and number_per_point is provided')
                number = number_per_point * len(lon)
            if number is not None:
                if number % len(lon) == 0:
                    number_per_point = int(number / len(lon))
                    if number_per_point > 1:
                        logger.debug(f'Seeding {number_per_point} elements per point')
                        lon = np.repeat(lon, number_per_point)
                        lat = np.repeat(lat, number_per_point)
                else:
                    raise ValueError(
                        'Lon and lat have length %s, but number is %s, which is not a multiple' %
                        (len(lon), number))
            number = len(lon)
        else:
            if number is None:
                if len(time) > 2:
                    number = len(time)  # Interpreting as time series
                else:
                    number = self.get_config('seed:number')
            lon = lon * np.ones(number)
            lat = lat * np.ones(number)

        if len(time) != number and len(time) > 1:
            if len(time) == 2:  # start -> end
                td = (time[1] - time[0]) / (number - 1
                                            )  # timestep between points
                time = [time[0] + i * td for i in range(number)]
            else:
                raise ValueError(
                    'Time array has length %s, must be 1, 2 or %s' %
                    (len(time), number))

        # Add radius / perturbation
        if radius.max() > 0:
            geod = pyproj.Geod(ellps='WGS84')
            if radius_type == 'gaussian':
                x = np.random.randn(np.sum(number)) * radius
                y = np.random.randn(np.sum(number)) * radius
                az = np.degrees(np.arctan2(x, y))
                dist = np.sqrt(x * x + y * y)
            elif radius_type == 'uniform':
                az = np.random.rand(np.sum(number)) * 360
                dist = np.sqrt(np.random.uniform(0, 1,
                                                 np.sum(number))) * radius
            lon, lat, az = geod.fwd(lon, lat, az, dist, radians=False)

        # If z is 'seafloor'
        if not 'z' in kwargs or kwargs['z'] is None:
            if 'seed:seafloor' in self._config:
                if self.get_config('seed:seafloor') is True:
                    kwargs['z'] = 'seafloor'
                    logger.debug('Seafloor is selected, neglecting z')
        if 'z' in kwargs and isinstance(kwargs['z'], str) \
                and kwargs['z'][0:8] == 'seafloor':
            # We need to fetch seafloor depth from reader
            seafloor_constant = self.get_config(
                'environment:constant:sea_floor_depth_below_sea_level')
            seafloor_fallback = self.get_config(
                'environment:fallback:sea_floor_depth_below_sea_level')
            if seafloor_constant is not None:
                env = {
                    'sea_floor_depth_below_sea_level':
                    np.array(seafloor_constant)
                }
            elif ('sea_floor_depth_below_sea_level'
                  in self.env.priority_list) or len(self.env._lazy_readers()):
                if not hasattr(self, 'time'):
                    self.time = time[0]
                env, env_profiles, missing = \
                    self.env.get_environment(['sea_floor_depth_below_sea_level'],
                                         time=time[0], lon=lon, lat=lat, z=0*lon)
            elif seafloor_fallback is not None:
                env = {
                    'sea_floor_depth_below_sea_level':
                    np.array(seafloor_fallback)
                }
            else:
                raise ValueError('A reader providing the variable '
                                 'sea_floor_depth_below_sea_level must be '
                                 'added before seeding elements at seafloor.')
            # Add M meters if given as 'seafloor+M'
            if len(kwargs['z']) > 8 and kwargs['z'][8] == '+':
                meters_above_seafloor = float(kwargs['z'][9::])
                logger.info('Seeding elements %f meters above seafloor' %
                            meters_above_seafloor)
            else:
                meters_above_seafloor = 0
            kwargs['z'] = \
                -env['sea_floor_depth_below_sea_level'].astype('float32') + meters_above_seafloor

        # Creating and scheduling elements
        elements = self.ElementType(lon=lon, lat=lat, **kwargs)
        time_array = np.array(time)
        self.schedule_elements(elements, time)

    @require_mode(mode=Mode.Ready)
    def seed_cone(self, lon, lat, time, radius=0, number=None, **kwargs):
        """Seed elements along a transect/cone between two points/times

        Arguments:
            lon: scalar or list with 2 elements [lon0, lon1]

            lat: scalar or list with 2 elements [lat0, lat]

            time: datetime or list with 2 elements [t0, t1]

            radius: scalar or list with 2 elements [r0, r1] Unit: meters

            number (int): The number of elements. If this is None, the number of
            elements is taken from configuration.

        Elements are seeded along a transect from
            (lon0, lat0) with uncertainty radius r0 at time t0, towards
            (lon1, lat1) with uncertainty radius r1 at time t1.
            If r0 != r1, the unceetainty radius is linearly changed along
            the transect, thus outlining a "cone".
        """

        if number is None:
            number = self.get_config('seed:number')
        if number == 1:
            raise ValueError(
                'For a cone, the number of elements must be at least 2 or more, given is 1'
            )

        lon = np.atleast_1d(lon).ravel()
        lat = np.atleast_1d(lat).ravel()
        radius = np.atleast_1d(radius).ravel()
        if len(lon) != len(lat):
            raise ValueError('Lon and lat must have same length (1 or 2)')
        elif len(lon) > 2:
            raise ValueError(
                'Lon and lat must have length 1 or 2, given length is %s' %
                (len(lon)))
        elif len(lon) == 1:
            lon = lon * np.ones(number)
            lat = lat * np.ones(number)
        elif len(lon) == 2:  # Segment from lon0,lat1 to lon1,lat2
            geod = pyproj.Geod(ellps='WGS84')
            lonin = lon
            latin = lat
            # Note that npts places points in-between start and end, and does not include these
            conelonlats = geod.npts(lon[0],
                                    lat[0],
                                    lon[1],
                                    lat[1],
                                    number,
                                    radians=False)
            lon, lat = zip(*conelonlats)

        if len(radius) > 2:
            raise ValueError('Seed radius must have length 1 or 2')
        elif len(radius) == 2:  # Linear increase from r0 to r1
            radius = np.linspace(radius[0], radius[1], number)

        if isinstance(time, list) and len(time) == 1:
            time = time[0]

        if hasattr(time, '__len__'):
            timespan = [time[0], time[-1]]
        else:
            timespan = [time, time]

        radius = radius.astype(np.float32)
        lonin = lonin if 'lonin' in locals() else [lon.min(), lon.max()]
        latin = latin if 'latin' in locals() else [lat.min(), lat.max()]

        self.seed_cone_arguments = {
            'lon': lonin,
            'lat': latin,
            'radius': [float(radius[0]), float(radius[-1])],
            'time': timespan,
            'number': number
        }

        # Make GeoJson seeding dict to be saved in netCDF metadata
        geo = geojson.LineString([(float(lonin[0]), float(latin[0])),
                                  (float(lonin[1]), float(latin[1]))])
        seed_defaults = self.get_configspec('seed')
        default_seed = {
            k.split(':')[-1]: seed_defaults[k]['value']
            for k in seed_defaults
        }
        if 'seafloor' in default_seed and default_seed['seafloor'] is True:
            default_seed['z'] = 'seafloor'
        default_seed = {
            **default_seed,
            **kwargs
        }  # Overwrite with explicitly provided values
        properties = {
            **default_seed, 'time': [str(timespan[0]),
                                     str(timespan[1])],
            'radius': [float(radius[0]), float(radius[-1])],
            'number': number
        }
        # convert array to string in case of array input to seed cone
        for key in properties.keys():
            if isinstance(properties[key], np.ndarray):
                properties[key] = np.array2string(properties[key])

        f = geojson.Feature(geometry=geo, properties=properties)
        self.seed_geojson.append(f)

        # Forwarding calculated cone points/radii to seed_elements
        self.seed_elements(lon=lon,
                           lat=lat,
                           time=time,
                           radius=radius,
                           number=number,
                           **kwargs)

    @require_mode(mode=Mode.Ready)
    def seed_from_geojson(self, gjson):
        """Under development"""
        try:
            gj = geojson.loads(gjson)
        except:
            raise ValueError('Could not load GeoJSON string: %s' % gjson)
        if not gj.is_valid:
            raise ValueError('GeoJSON string is not valid: %s' % gj.errors())
        # Assuming temporally that g is a Feature, and not a FeatureCollection
        properties = gj['properties']
        if 'time' not in properties:
            raise ValueError('Property "time" is not available')
        kwargs = {}
        for prop in properties:
            if prop == 'time':
                t = properties['time']
                if isinstance(t, list):
                    time = [
                        datetime.fromisoformat(t[0].replace("Z", "+00:00")),
                        datetime.fromisoformat(t[1].replace("Z", "+00:00"))
                    ]
                    time = [t.replace(tzinfo=None) for t in time]
                else:
                    time = datetime.fromisoformat(t.replace("Z", "+00:00"))
                    time = time.replace(tzinfo=None)
            else:
                kwargs[prop] = properties[prop]

        geometry = gj['geometry']

        if geometry['type'] == 'Polygon':
            coords = list(geojson.utils.coords(gj))
            lon, lat = zip(*[(c[0], c[1]) for c in coords])
            self.seed_within_polygon(lons=lon, lats=lat, time=time, **kwargs)
        elif geometry['type'] == 'LineString':
            coords = list(geojson.utils.coords(gj))
            lon, lat = zip(*[(c[0], c[1]) for c in coords])
            self.seed_cone(lon=lon, lat=lat, time=time, **kwargs)
        elif geometry['type'] == 'Point':
            coords = list(geojson.utils.coords(gj))
            lon, lat = zip(*[(c[0], c[1]) for c in coords])
            self.seed_elements(lon=lon, lat=lat, time=time, **kwargs)
        else:
            raise ValueError('Not yet implemented')

    @require_mode(mode=Mode.Ready)
    def seed_repeated_segment(self,
                              lons,
                              lats,
                              start_time,
                              end_time,
                              time_interval=None,
                              number_per_segment=None,
                              total_number=None,
                              **kwargs):
        """Seed elements repeatedly in time along a segment.

        The segment goes from lon[0],lat[0] to lon[1],lat[1].

        The number of elements should be proved as either:

        1) number_per_segment, in which case total number of elements is number_per_segment * len(times), or

        2) total_number, in which case the number of elements per segment is: total_number / len(times).
           Any extra elements are duplicated along at the first segment.

        """

        numtimes = int((end_time - start_time).total_seconds() /
                       time_interval.total_seconds() + 1)
        times = [start_time + i * time_interval for i in range(numtimes)]

        geod = pyproj.Geod(ellps='WGS84')
        if number_per_segment is None:
            number_per_segment = int(np.floor(total_number / numtimes))

        s_lonlats = geod.npts(lons[0],
                              lats[0],
                              lons[1],
                              lats[1],
                              number_per_segment,
                              radians=False)
        slon, slat = list(zip(*s_lonlats))
        slon = np.atleast_1d(slon)
        slat = np.atleast_1d(slat)

        lon, time = np.meshgrid(slon, times)
        lat, time = np.meshgrid(slat, times)
        lon = lon.ravel()
        lat = lat.ravel()
        time = time.ravel()

        if total_number is not None:
            additional_elements = total_number - len(lon.ravel())
            logger.info('Repeating the %d last points, to obtain %d elements' %
                  (additional_elements, total_number))
            lon = np.concatenate((lon, lon[-additional_elements::]))
            lat = np.concatenate((lat, lat[-additional_elements::]))
            time = np.concatenate((time, time[-additional_elements::]))

        self.seed_elements(lon=lon, lat=lat, time=time, **kwargs)

    @require_mode(mode=Mode.Ready)
    def seed_within_polygon(self, lons, lats, number=None, **kwargs):
        """Seed a number of elements within given polygon.

        Arguments:
            lon: array of longitudes

            lat: array of latitudes

            number: int, number of elements to be seeded

            kwargs: keyword arguments containing properties/attributes and
            values corresponding to the actual particle type (ElementType).
            These are forwarded to method seed_elements(). All properties
            for which there are no default value must be specified.

        """
        if number == 0:
            return

        if number is None:
            number = self.get_config('seed:number')

        lons = np.asarray(lons)
        lats = np.asarray(lats)
        if len(lons) < 3:
            logger.info('At least three points needed to make a polygon')
            return
        if len(lons) != len(lats):
            raise ValueError('lon and lat arrays must have same length.')
        poly = Polygon(list(zip(lons, lats)), closed=True)
        # Place N points within the polygons
        proj = pyproj.Proj('+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f '
                           '+lon_0=%f +R=6370997.0 +units=m +ellps=WGS84' %
                           (lats.min(), lats.max(),
                            (lats.min() + lats.max()) / 2,
                            (lons.min() + lons.max()) / 2))
        lonlat = poly.get_xy()
        lon = lonlat[:, 0]
        lat = lonlat[:, 1]
        x, y = proj(lon, lat)
        area = 0.0
        for i in range(-1, len(x) - 1):
            area += x[i] * (y[i + 1] - y[i - 1])
        area = abs(area) / 2

        # Make points, evenly distributed
        deltax = np.sqrt(area / number)
        lonpoints = np.array([])
        latpoints = np.array([])
        lonlat = poly.get_xy()
        lon = lonlat[:, 0]
        lat = lonlat[:, 1]
        x, y = proj(lon, lat)
        xvec = np.linspace(x.min() + deltax / 2,
                           x.max() - deltax / 2,
                           int((x.max() - x.min()) / deltax))
        yvec = np.linspace(y.min() + deltax / 2,
                           y.max() - deltax / 2,
                           int((y.max() - y.min()) / deltax))
        x, y = np.meshgrid(xvec, yvec)
        lon, lat = proj(x, y, inverse=True)
        lon = lon.ravel()
        lat = lat.ravel()
        points = np.c_[lon, lat]
        ind = Path(poly.xy).contains_points(points)
        if not any(ind):  # No elements are inside, we seed on border
            lonpoints = np.append(lonpoints, lons[0:number])
            latpoints = np.append(latpoints, lats[0:number])
        else:
            lonpoints = np.append(lonpoints, lon[ind])
            latpoints = np.append(latpoints, lat[ind])
        if len(ind) == 0:
            logger.info('Small or irregular polygon, using center point.')
            lonpoints = np.atleast_1d(np.mean(lons))
            latpoints = np.atleast_1d(np.mean(lats))
        # Truncate if too many
        # NB: should also repeat some points, if too few
        lonpoints = lonpoints[0:number]
        latpoints = latpoints[0:number]
        while True:
            if len(lonpoints) < number:
                # If number of positions is smaller than requested,
                # we duplicate the first ones
                missing = number - len(lonpoints)
                lonpoints = np.append(lonpoints, lonpoints[0:missing])
                latpoints = np.append(latpoints, latpoints[0:missing])
            else:
                break

        # Finally seed at calculated positions
        self.seed_elements(lonpoints, latpoints, number=number, **kwargs)

    @require_mode(mode=Mode.Ready)
    def seed_from_wkt(self, wkts, time, **kwargs):
        """Seeds elements within (multi)polygons from WKT"""
        from shapely import wkt
        import geopandas as gpd  # Local import as this is not a requirement
        df = pd.DataFrame({'id': [0], 'geometry': wkts})
        df['geometry'] = df['geometry'].apply(wkt.loads)
        gdf = gpd.GeoDataFrame(df, geometry='geometry')
        if gdf.crs is None:
            gdf = gdf.set_crs(4326, allow_override=True)
        self.seed_from_geopandas(gdf, time, **kwargs)

    @require_mode(mode=Mode.Ready)
    def seed_from_shapefile(self,
                            shapefile,
                            time,
                            **kwargs):
        """Seeds elements within polygons of a shapefile"""

        import geopandas as gpd  # Local import as this is not a requirement
        geodataframe = gpd.read_file(shapefile)
        self.seed_from_geopandas(geodataframe, time, **kwargs)

    @require_mode(mode=Mode.Ready)
    def seed_from_geopandas(self,
                            geodataframe,
                            time,
                            **kwargs):
        """Seeds elements within polygons of a GeoPandas DataFrame"""

        g = geodataframe

        geometry_types = g.geom_type.unique()
        if geometry_types == ['Point']:
            logger.debug('Geodataframe contains only points, calling seed_elements()')
            self.seed_elements(lon=g.geometry.x, lat=g.geometry.y, time=time, **kwargs)
            return

        # Below only Polygon and Multiplolygon is assumed.
        # TODO:
        # - handling of linestring
        # - raise error if mixture of point, linestring and polygons

        if 'number' not in kwargs:
            number = self.get_config('seed:number')
        else:
            number = kwargs.pop('number')

        ga = g.to_crs({'proj':'cea'})  # Equal area projection for area calculation
        g_lonlat = g.to_crs({'proj':'lonlat'})  # Lonlat projection to get lon and lat
        
        areas = np.array([p.area for p in ga.geometry.explode(index_parts=False)])
        logger.info(f'Seeding {number} elements within {len(areas)} polygons')

        total_area = np.sum(areas)
        number_per_polygon = np.array([round(a/total_area*number) for a in areas])
        # If numbers do not sum to total number, we add/remove 1 particles for the first polygons
        remaining = int(number - sum(number_per_polygon))
        if abs(remaining) <= len(number_per_polygon):
            number_per_polygon[0:abs(remaining)] += np.sign(remaining)
        elif abs(remaining) > len(number_per_polygon):
            raise ValueError('Should not happen')
        
        for e, (n, polygon) in enumerate(zip(number_per_polygon, g_lonlat.geometry.explode(index_parts=False))):
            logger.info(f'Seeding {n} elements within polygon number {e+1} of area {areas[e]/1e6} km2')
            lons, lats = polygon.exterior.coords.xy
            self.seed_within_polygon(lons=lons,
                                     lats=lats,
                                     number=n,
                                     time=time,
                                     **kwargs)

    # @require_mode(mode=Mode.Ready)
    # def seed_from_shapefile_old(self,
    #                         shapefile,
    #                         number,
    #                         layername=None,
    #                         featurenum=None,
    #                         **kwargs):
    #     """Seeds elements within contours read from a shapefile
        
    #         Obsolete, as new method based on geopandas is simpler"""

    #     try:
    #         from osgeo import ogr, osr
    #     except Exception as e:
    #         logger.warning(e)
    #         raise ValueError('OGR library is needed to read shapefiles.')

    #     if 'timeformat' in kwargs:
    #         # Recondstructing time from filename, where 'timeformat'
    #         # is forwarded to datetime.strptime()
    #         kwargs['time'] = datetime.strptime(os.path.basename(shapefile),
    #                                            kwargs['timeformat'])
    #         del kwargs['timeformat']

    #     num_seeded_before = self.num_elements_scheduled()

    #     targetSRS = osr.SpatialReference()
    #     targetSRS.ImportFromEPSG(4326)
    #     try:
    #         s = ogr.Open(shapefile)
    #     except:
    #         s = shapefile

    #     for layer in s:
    #         if layername is not None and layer.GetName() != layername:
    #             logger.info('Skipping layer: ' + layer.GetName())
    #             continue
    #         else:
    #             logger.info('Seeding for layer: %s (%s features)' %
    #                         (layer.GetDescription(), layer.GetFeatureCount()))

    #         coordTrans = osr.CoordinateTransformation(layer.GetSpatialRef(),
    #                                                   targetSRS)

    #         if featurenum is None:
    #             featurenum = range(1, layer.GetFeatureCount() + 1)
    #         else:
    #             featurenum = np.atleast_1d(featurenum)
    #         if max(featurenum) > layer.GetFeatureCount():
    #             raise ValueError('Only %s features in layer.' %
    #                              layer.GetFeatureCount())

    #         # Loop first through all features to determine total area
    #         layer.ResetReading()
    #         area_srs = osr.SpatialReference()
    #         area_srs.ImportFromEPSG(3857)
    #         areaTransform = osr.CoordinateTransformation(
    #             layer.GetSpatialRef(), area_srs)

    #         areas = np.zeros(len(featurenum))
    #         for i, f in enumerate(featurenum):
    #             feature = layer.GetFeature(f - 1)  # Note 1-indexing, not 0
    #             if feature is not None:
    #                 gom = feature.GetGeometryRef().Clone()
    #                 gom.Transform(areaTransform)
    #                 areas[i] = gom.GetArea()

    #         total_area = np.sum(areas)
    #         layer.ResetReading()  # Rewind to first layer
    #         logger.info('Total area of all polygons: %s m2' % total_area)
    #         # Find number of points per polygon
    #         numbers = np.round(number * areas / total_area).astype(int)
    #         numbers[numbers.argmax()] += int(number - sum(numbers))

    #         for i, f in enumerate(featurenum):
    #             feature = layer.GetFeature(f - 1)
    #             if feature is None:
    #                 continue
    #             num_elements = numbers[i]
    #             geom = feature.GetGeometryRef()
    #             logger.info(f'\tSeeding {num_elements} elements within polygon number {featurenum[i]} of area {areas[i]} m3')
    #             try:
    #                 geom.Transform(coordTrans)
    #             except Exception as e:
    #                 logger.warning('Could not transform coordinates:')
    #                 logger.warning(e)
    #                 pass
    #             #b = geom.GetBoundary()
    #             #if b is not None:
    #             #    points = b.GetPoints()
    #             #    lons = [p[0] for p in points]
    #             #    lats = [p[1] for p in points]
    #             #else:
    #             # Alternative if OGR is not built with GEOS support
    #             r = geom.GetGeometryRef(0)
    #             lons = [r.GetY(j) for j in range(r.GetPointCount())]
    #             lats = [r.GetX(j) for j in range(r.GetPointCount())]

    #             self.seed_within_polygon(lons=lons,
    #                                      lats=lats,
    #                                      number=num_elements,
    #                                      **kwargs)

    @require_mode(mode=Mode.Ready)
    def seed_letters(self, text, lon, lat, time, number, scale=1.2):
        """Seed elements within text polygons"""
        from matplotlib.font_manager import FontProperties
        fp = FontProperties(family='Bitstream Vera Sans', weight='bold')
        pol = matplotlib.textpath.TextPath((lon, lat),
                                           text,
                                           size=1 * scale,
                                           prop=fp)
        patch = matplotlib.patches.PathPatch(pol,
                                             facecolor='none',
                                             edgecolor='black',
                                             transform=ccrs.PlateCarree())
        po = patch.get_path().to_polygons()
        for p in po:
            self.seed_within_polygon(lons=p[:, 0],
                                     lats=p[:, 1],
                                     number=number,
                                     time=time)

    @require_mode(mode=Mode.Ready)
    def seed_from_ladim(self, ladimfile, roms):
        """Seed elements from ladim \\*.rls text file: [time, x, y, z, name]"""

        data = np.loadtxt(ladimfile,
                          dtype={
                              'names': ('time', 'x', 'y', 'z'),
                              'formats': ('S20', 'f4', 'f4', 'f4')
                          },
                          usecols=(0, 1, 2, 3))

        time = [datetime.strptime(t, "%Y-%m-%dT%H") for t in data['time']]
        time = np.array(time)

        lon, lat = roms.xy2lonlat(data['x'], data['y'])
        z = -data['z']

        logger.info('Seeding %i elements from %s:' % (len(lon), ladimfile))
        logger.info('    Lons: %f to %f' % (lon.min(), lon.max()))
        logger.info('    Lats: %f to %f' % (lat.min(), lat.max()))
        logger.info('    Depths: %f to %f' % (z.min(), z.max()))
        logger.info('    Time: %s to %s' % (time.min(), time.max()))
        elements = self.ElementType(lon=lon, lat=lat, z=-z)

        self.schedule_elements(elements, time)

    def horizontal_diffusion(self):
        """Move elements with random walk according to given horizontal diffuivity."""
        D = self.get_config('drift:horizontal_diffusivity')
        if D == 0:
            logger.debug('Horizontal diffusivity is 0, no random walk.')
            return
        if self.num_elements_active() == 0:
            logger.debug('No active elements, skipping horizontal diffusivity.')
            return
        dt = np.abs(self.time_step.total_seconds())
        x_vel = self.elements.moving * np.sqrt(2 * D / dt) * np.random.normal(
            scale=1, size=self.num_elements_active())
        y_vel = self.elements.moving * np.sqrt(2 * D / dt) * np.random.normal(
            scale=1, size=self.num_elements_active())
        speed = np.sqrt(x_vel * x_vel + y_vel * y_vel)
        logger.debug(
            'Moving elements according to horizontal diffusivity of %s, with speeds between %s and %s m/s'
            % (D, speed.min(), speed.max()))
        self.update_positions(x_vel, y_vel)

    def deactivate_elements(self, indices, reason='deactivated'):
        """Schedule deactivated particles for deletion (at end of step)"""
        if any(indices) is False:
            return
        if reason not in self.status_categories:
            self.status_categories.append(reason)
            logger.debug('Added status %s' % (reason))
        reason_number = self.status_categories.index(reason)
        #if not hasattr(self.elements.status, "__len__"):
        if len(np.atleast_1d(self.elements.status)) == 1:
            status = self.elements.status.item()
            self.elements.status = np.zeros(self.num_elements_active())
            self.elements.status.fill(status)
        # Deactivate elements, if they have not already been deactivated
        self.elements.status[indices & (self.elements.status ==0)] = \
            reason_number
        self.elements.moving[indices] = 0
        logger.debug('%s elements scheduled for deactivation (%s)' %
                     (np.sum(indices), reason))
        logger.debug(
            '\t(z: %f to %f)' %
            (self.elements.z[indices].min(), self.elements.z[indices].max()))

    def remove_deactivated_elements(self):
        """Moving deactivated elements from self.elements
        to self.elements_deactivated."""

        # All particles scheduled for deletion
        indices = (self.elements.status != 0)
        #try:
        #    len(indices)
        #except:
        if len(indices) == 0 or np.sum(indices) == 0:
            logger.debug('No elements to deactivate')
            return  # No elements scheduled for deactivation
        # Basic, but some more housekeeping will be required later
        self.elements.move_elements(self.elements_deactivated, indices)
        logger.debug('Removed %i elements.' % (np.sum(indices)))
        if hasattr(self, 'environment'):
            self.environment = self.environment[~indices]
            logger.debug('Removed %i values from environment.' %
                         (np.sum(indices)))
        if hasattr(self, 'environment_profiles') and \
                self.environment_profiles is not None:
            for varname, profiles in self.environment_profiles.items():
                logger.debug('remove items from profile for ' + varname)
                if varname != 'z':
                    self.environment_profiles[varname] = \
                        profiles[:, ~indices]
            logger.debug('Removed %i values from environment_profiles.' %
                         (np.sum(indices)))
            #if self.num_elements_active() == 0:
            #    raise ValueError('No more active elements.')  # End simulation

    @require_mode(mode=Mode.Run, post_next_mode=True)
    def run(self,
            time_step=None,
            steps=None,
            time_step_output=None,
            duration=None,
            end_time=None,
            outfile=None,
            export_variables=None,
            export_buffer_length=100,
            stop_on_error=False):
        """Start a trajectory simulation, after initial configuration.

        Performs the main loop:
            - Obtain environment data for positions of all particles.
            - Call method 'update' to update (incl advect) particle properties.
        until one of the following conditions are met:
            - Maximum number of steps are reached
            - A needed variable can not be obtained by any reader
                (outside spatial/temporal domain) and has no fallback
                (default) value.
            - All particles have been deactivated (e.g. by stranding)
            - Occurance of any error, whose trace will be output to terminal.

        Before starting a model run, readers must be added for all
        required variables, unless fallback values have been specified.
        Some particles/elements must have been scheduled for seeding, and the
        run will start at the time when the first element has been scheduled..

        Arguments:
            time_step: interval between particles updates, in seconds or as
                timedelta. Default: 3600 seconds (1 hour)
            time_step_output: Time step at which element properties are stored
                and eventually written to file.
                Timedelta object or seconds.
                Default: same as time_step, meaning that all steps are stored
            The length of the simulation is specified by defining one
                (and only one) of the following parameters:
                - steps: integer, maximum number of steps. End of simulation
                will be self.start_time + steps*self.time_step
                - duration: timedelta defining the length of the simulation
                - end_time: datetime object defining the end of the simulation
            export_variables: list of variables and parameter names to be
                saved to file. Default is None (all variables are saved)
        """

        # Exporting software and hardware specification, for possible debugging
        logger.debug(opendrift.versions())

        self.timer_end('configuration')
        self.timer_start('preparing main loop')

        if self.num_elements_scheduled() == 0:
            raise ValueError('Please seed elements before starting a run.')
        self.elements = self.ElementType()

        # Export seed_geojson as FeatureCollection string
        self.add_metadata('seed_geojson',
                          geojson.FeatureCollection(self.seed_geojson))

        if outfile is None and export_buffer_length is not None:
            logger.debug('No output file is specified, '
                         'neglecting export_buffer_length')
            export_buffer_length = None

        # Some cleanup needed if starting from imported state
        if self.steps_calculation >= 1:
            self.steps_calculation = 0

        ########################
        # Simulation time step
        ########################
        if time_step is None:
            time_step = timedelta(
                minutes=self.get_config('general:time_step_minutes'))
        if type(time_step) is not timedelta:
            # Time step may be given in seconds, as alternative to timedelta
            time_step = timedelta(seconds=time_step)
        self.time_step = time_step
        if time_step_output is None:
            time_step_output = self.get_config(
                'general:time_step_output_minutes')
            if time_step_output is None:
                self.time_step_output = self.time_step
            else:
                self.time_step_output = timedelta(minutes=time_step_output)
        else:
            if type(time_step_output) is timedelta:
                self.time_step_output = time_step_output
            else:
                self.time_step_output = timedelta(seconds=time_step_output)
            if self.time_step_output.days >= 0 and self.time_step.days < 0:
                self.time_step_output = -self.time_step_output

        time_step_ratio = self.time_step_output.total_seconds() / \
            self.time_step.total_seconds()
        if time_step_ratio < 1:
            raise ValueError('Output time step must be equal or larger '
                             'than calculation time step.')
        if not time_step_ratio.is_integer():
            raise ValueError('Ratio of calculation and output time steps '
                             'must be an integer - given ratio is %s' %
                             time_step_ratio)
        ########################
        # Simulation duration
        ########################
        if time_step.days < 0:
            logger.info(
                'Backwards simulation, starting from last seeded element')
            self.start_time = self.elements_scheduled_time.max()
        if (duration is not None and end_time is not None) or \
            (duration is not None and steps is not None) or \
                (steps is not None and end_time is not None):
            raise ValueError('Only one of "steps", "duration" and "end_time" '
                             'may be provided simultaneously')
        if duration is None and end_time is None:
            if steps is not None:
                duration = steps * self.time_step
            else:
                for reader in self.env.readers.values():
                    if reader.end_time is not None:
                        if end_time is None:
                            end_time = reader.end_time
                        else:
                            end_time = min(end_time, reader.end_time)
                    logger.info('Duration, steps or end time not specified, '
                                'running until end of first reader: %s' %
                                (end_time))
        if duration is None:
            duration = end_time - self.start_time

        if time_step.days < 0 and duration.days >= 0:
            # Duration shall also be negative for backwards run
            duration = -duration

        if np.sign(duration.total_seconds()) * np.sign(
                time_step.total_seconds()) < 0:
            raise ValueError(
                "Time step must be negative if duration is negative.")

        ratio_duration_output = duration/self.time_step_output
        if not ratio_duration_output.is_integer():
            initial_duration = duration
            duration = np.ceil(ratio_duration_output)*self.time_step_output
            logger.warning('Simulation end is not at an output time step. '
                           f'Extending duration from {initial_duration} to {duration}')

        self.expected_steps_output = duration.total_seconds() / \
            self.time_step_output.total_seconds() + 1  # Includes start and end
        self.expected_steps_calculation = duration.total_seconds() / \
            self.time_step.total_seconds()
        self.expected_steps_output = int(self.expected_steps_output)
        self.expected_steps_calculation = int(self.expected_steps_calculation)
        self.expected_end_time = self.start_time + self.expected_steps_calculation * self.time_step

        ##############################################################
        # Prepare readers for the requested simulation domain/time
        ##############################################################
        max_distance = \
                self.get_config('drift:max_speed')*self.expected_steps_calculation * \
            np.abs(self.time_step.total_seconds())
        deltalat = max_distance / 111000.
        deltalon = deltalat / np.cos(
            np.radians(np.mean(self.elements_scheduled.lat)))
        # TODO: extent should ideally be a general polygon, not only lon/lat-min/max
        # TODO: Should also take into account eventual lifetime of elements
        simulation_extent = [
            np.maximum(-360,
                       self.elements_scheduled.lon.min() - deltalon),
            np.maximum(-89,
                       self.elements_scheduled.lat.min() - deltalat),
            np.minimum(360,
                       self.elements_scheduled.lon.max() + deltalon),
            np.minimum(89,
                       self.elements_scheduled.lat.max() + deltalat)
        ]
        if simulation_extent[2] == 360 and simulation_extent[0] < 0:
            simulation_extent[0] = 0

        logger.debug(
            'Finalizing environment and preparing readers for simulation coverage (%s) and time (%s to %s)'
            % (simulation_extent, self.start_time, self.expected_end_time))

        # Store expected simulation extent, to check if new readers have coverage
        self.simulation_extent = simulation_extent
        self.env.finalize(self.simulation_extent)

        ####################################################################
        # Preparing history array for storage in memory and eventually file
        ####################################################################
        if export_buffer_length is None:
            self.export_buffer_length = self.expected_steps_output
        else:
            self.export_buffer_length = np.minimum(export_buffer_length, self.expected_steps_output)

        if self.time_step.days < 0:
            # For backwards simulation, we start at last seeded element
            logger.info('Backwards simulation, starting at '
                        'time of last seeded element')
            self.time = self.elements_scheduled_time.max()
            # Flipping ID array, so that lowest IDs are released first
            self.elements_scheduled.ID = \
                np.flipud(self.elements_scheduled.ID)
        else:
            # Forward simulation, start time has been set when seeding
            self.time = self.start_time

        # Add the output variables which are always required
        if export_variables is not None:
            export_variables = list(
                set(export_variables + ['lon', 'lat', 'status']))
        self.export_variables = export_variables

        # Create Xarray Dataset to hold result
        coords = {  # Initialize for the part fitting in memory
            'trajectory': ('trajectory', np.arange(len(self.elements_scheduled)),
                {'cf_role': 'trajectory_id', 'dtype': np.int32}),
            'time': ('time', pd.date_range(self.start_time, periods=self.export_buffer_length,
                                           freq=self.time_step_output),
                     {'standard_name': 'time', 'long_name': 'time', 'dtype': np.float64})}
        shape = (len(coords['trajectory'][1]), len(coords['time'][1]))
        dims = ('trajectory', 'time')  # Presently, but shall also allow single dimension

        element_vars = {}
        default_dtype = np.float32  # Allows NaN (in contrast to np.int)
        for varname, attrs in self.ElementType.variables.items():
            if varname == 'ID' or (self.export_variables is not None and
                                   varname not in self.export_variables):
                continue
            attrs = {k:v for k,v in attrs.copy().items() if k not in ['seed', 'default']}
            element_vars[varname] = (dims,
                                     np.full(shape=shape, fill_value=np.nan, dtype=default_dtype),
                                     attrs)

        environment_vars = {varname: (dims, np.full(shape=shape, fill_value=np.nan, dtype=default_dtype))
                            for varname,var in self.required_variables.items()
                            if self.export_variables is None or varname in self.export_variables}
        global_attributes = {
            'Conventions': 'CF-1.11, ACDD-1.3',
            'standard_name_vocabulary': 'CF Standard Name Table v85',
            'featureType': 'trajectory',
            'title': 'OpenDrift trajectory simulation',
            'summary': 'Output from simulation with OpenDrift framework',
            'keywords': 'trajectory, drift, lagrangian, simulation',
            'history': 'Created ' + str(datetime.now()),
            'date_created': datetime.now().isoformat(),
            'source': 'Output from simulation with OpenDrift',
            'model_url': 'https://github.com/OpenDrift/opendrift',
            'opendrift_class': self.__class__.__name__,
            'opendrift_module': self.__class__.__module__,
            'readers': str(self.env.readers.keys()),
            'time_coverage_start': str(self.start_time),
            'time_step_calculation': str(self.time_step),
            'time_step_output': str(self.time_step_output),
            }

        # Add config settings as global attributes
        for key in self._config:
            value = self.get_config(key)
            if isinstance(value, (bool, type(None))):
                value = str(value)
            global_attributes['config_' + key] = value

        # TODO: add/remove metadata-dict

        self.result = xr.Dataset(coords=coords, data_vars=element_vars | environment_vars, attrs=global_attributes)
        logger.debug(f'Initial self.result, size {self.result.sizes}')

        if self.origin_marker is not None and 'origin_marker' in self.result.data_vars:
            self.result['origin_marker'] = self.result.origin_marker.assign_attrs(
                {'flag_values': np.arange(len(self.origin_marker)).astype(self.ElementType.variables['origin_marker']['dtype']),
                 'flag_meanings': " ".join(self.origin_marker.values())
                })

        if outfile is not None:
            self.io_init(outfile)
            self.outfile = outfile
        else:
            self.outfile = None

        # Move point seeded on land to ocean
        if self.get_config('seed:ocean_only') is True and \
            ('land_binary_mask' in self.required_variables):
            #('land_binary_mask' not in self.fallback_values) and \
            self.timer_start('preparing main loop:moving elements to ocean')
            self.elements_scheduled.lon, self.elements_scheduled.lat = \
                self.closest_ocean_points(self.elements_scheduled.lon,
                                          self.elements_scheduled.lat)
            self.timer_end('preparing main loop:moving elements to ocean')

        #############################
        # Check validity domain
        #############################
        validity_domain = [
            self.get_config('drift:deactivate_west_of'),
            self.get_config('drift:deactivate_east_of'),
            self.get_config('drift:deactivate_south_of'),
            self.get_config('drift:deactivate_north_of')
        ]
        if validity_domain == [None, None, None, None]:
            self.validity_domain = None
        else:
            self.validity_domain = validity_domain

        #############################
        # Model specific preparation
        #############################
        self.prepare_run()

        ##########################
        # Main loop
        ##########################
        self.add_metadata('simulation_time', datetime.now())
        self.timer_end('preparing main loop')
        self.timer_start('main loop')
        self.memory_usage = np.array([])
        for i in range(self.expected_steps_calculation):
            self.memory_usage = np.append(self.memory_usage, psutil.virtual_memory().used / (1024.0**3))
            try:
                # Release elements
                self.release_elements()

                if self.num_elements_active(
                ) == 0 and self.num_elements_scheduled() > 0:
                    logger.info(
                        'No active but %s scheduled elements, skipping timestep %s (%s)'
                        % (self.num_elements_scheduled(),
                           self.steps_calculation + 1 , self.time))
                    self.state_to_buffer()  # Append status to history array
                    self.steps_calculation += 1
                    if self.time is not None:
                        self.time = self.time + self.time_step
                    continue

                self.interact_with_seafloor()

                if self.show_continuous_performance is True:
                    logger.info(self.performance())
                # Display time to terminal
                logger.debug('===================================' * 2)
                logger.info('%s - step %i of %i - %i active elements '
                            '(%i deactivated)' %
                            (self.time, self.steps_calculation + 1,
                             self.expected_steps_calculation,
                             self.num_elements_active(),
                             self.num_elements_deactivated()))
                logger.debug('%s elements scheduled.' %
                             self.num_elements_scheduled())
                logger.debug('===================================' * 2)
                if len(self.elements.lon) > 0:
                    lonmin = self.elements.lon.min()
                    lonmax = self.elements.lon.max()
                    latmin = self.elements.lat.min()
                    latmax = self.elements.lat.max()
                    zmin = self.elements.z.min()
                    zmax = self.elements.z.max()
                    if latmin == latmax:
                        logger.debug('\t\tlatitude =  %s' % (latmin))
                    else:
                        logger.debug('\t\t%s <- latitude  -> %s' %
                                     (latmin, latmax))
                    if lonmin == lonmax:
                        logger.debug('\t\tlongitude = %s' % (lonmin))
                    else:
                        logger.debug('\t\t%s <- longitude -> %s' %
                                     (lonmin, lonmax))
                    if zmin == zmax:
                        logger.debug('\t\tz = %s' % (zmin))
                    else:
                        logger.debug('\t\t%s   <- z ->   %s' % (zmin, zmax))
                    logger.debug('---------------------------------')

                self.environment, self.environment_profiles, missing = \
                    self.env.get_environment(list(self.required_variables),
                                         self.time,
                                         self.elements.lon,
                                         self.elements.lat,
                                         self.elements.z,
                                         self.required_profiles,
                                         self.profiles_depth)

                self.store_previous_variables()

                self.calculate_missing_environment_variables()

                if any(missing):
                    self.report_missing_variables()

                self.interact_with_coastline()

                self.interact_with_seafloor()

                self.deactivate_elements(missing, reason='missing_data')

                self.state_to_buffer()  # Append status to history array

                self.increase_age_and_retire()

                self.remove_deactivated_elements()

                # Propagate one timestep forwards
                self.steps_calculation += 1

                if self.num_elements_active(
                ) == 0 and self.num_elements_scheduled() == 0:
                    raise ValueError(
                        'No more active or scheduled elements, quitting.')

                # Store location, in case elements shall be moved back
                self.store_present_positions()

                #####################################################
                if self.num_elements_active() > 0:
                    logger.debug('Calling %s.update()' % type(self).__name__)
                    self.timer_start('main loop:updating elements')
                    self.update()
                    self.timer_end('main loop:updating elements')
                else:
                    logger.info('No active elements, skipping update() method')
                #####################################################

                self.horizontal_diffusion()

                if self.num_elements_active(
                ) == 0 and self.num_elements_scheduled() == 0:
                    raise ValueError(
                        'No active or scheduled elements, quitting simulation')

                logger.debug('%s active elements (%s deactivated)' %
                             (self.num_elements_active(),
                              self.num_elements_deactivated()))
                # Updating time
                if self.time is not None:
                    self.time = self.time + self.time_step

            except Exception as e:
                message = ('The simulation stopped before requested '
                           'end time was reached.')
                logger.warning(message)
                self.store_message(message)
                logger.info('========================')
                logger.info('End of simulation:')
                logger.info(e)
                logger.info(traceback.format_exc())
                logger.info(self.get_messages())
                if not hasattr(self, 'environment'):
                    sys.exit('Simulation aborted. ' + self.get_messages())
                logger.info('========================')
                if stop_on_error is True:
                    sys.exit('Stopping on error. ' + self.get_messages())
                if self.steps_calculation <= 1:
                    raise ValueError('Simulation stopped within '
                                     'first timestep. ' + self.get_messages())
                break

        self.timer_end('main loop')
        self.timer_start('cleaning up')
        logger.debug('Cleaning up')

        self.interact_with_coastline(final=True)
        self.timer_end('cleaning up')
        self.timer_end('total time')
        self.state_to_buffer(final=True)  # Append final status to buffer

        ## Add any other data to the result here.
        self.post_run()

        if outfile is not None:
            logger.debug('Finalising and closing output file: %s' % outfile)
            self.io_close()
        else:
            if self.num_elements_scheduled() > 0:
                logger.info(f'Removing {self.num_elements_scheduled()} unseeded elements')
                seeded_indices = [n for n in np.arange(self.num_elements_total())
                             if n not in self.elements_scheduled.ID]
                self.result = self.result.isel(trajectory=seeded_indices)

        # Remove any elements scheduled for deactivation during last step
        self.remove_deactivated_elements()

        if export_buffer_length is None:
            pass  # TODO - do this for self.result
        else:  # If output has been flushed to file during run, we
            # need to reimport from file to get all data in memory
            del self.environment
            if hasattr(self, 'environment_profiles'):
                del self.environment_profiles
            self.result = xr.open_dataset(outfile)

        return self.result

    def increase_age_and_retire(self):
        """Increase age of elements, and retire if older than config setting."""
        # Increase age of elements
        self.elements.age_seconds += self.time_step.total_seconds()

        # Deactivate elements that exceed a certain age
        if self.get_config('drift:max_age_seconds') is not None:
            self.deactivate_elements(
                self.elements.age_seconds
                >= self.get_config('drift:max_age_seconds'),
                reason='retired')

        # Deacticate any elements outside validity domain set by user
        if self.validity_domain is not None:
            W, E, S, N = self.validity_domain
            if W is not None:
                self.deactivate_elements(self.elements.lon < W,
                                         reason='outside')
            if E is not None:
                self.deactivate_elements(self.elements.lon > E,
                                         reason='outside')
            if S is not None:
                self.deactivate_elements(self.elements.lat < S,
                                         reason='outside')
            if N is not None:
                self.deactivate_elements(self.elements.lat > N,
                                         reason='outside')

    def state_to_buffer(self, final=False):

        if pd.to_datetime(self.time) in self.result.time:  # Output time step
            ID_ind = self.elements.ID
            element_ind = range(len(ID_ind))
            insert_time = pd.to_datetime(self.time)
        else:  # Deactivated elements must be written even if no output timestep
            deactivated = np.where(self.elements.status != 0)[0]
            if len(deactivated) == 0 and final is False:
                return  # No deactivated elements this sub-timestep
            ID_ind = self.elements.ID[deactivated]
            element_ind = deactivated
            insert_time = self.result.time.sel(time=self.time, method='backfill').values

        for var in self.result.data_vars:
            for source in (self.elements, self.environment):
                d = getattr(source, var, None)
                if d is not None:
                    self.result[var].loc[{'time': insert_time,
                                          'trajectory': ID_ind}] = d[element_ind]

        buffer_is_full = (pd.to_datetime(self.time) == self.result.time[-1]).values

        if final is True or buffer_is_full:
            logger.debug('Updating minval and maxval')
            # Update min and max values
            for varname, var in self.result.data_vars.items():
                if varname in ['status']:
                    continue
                minval = var.min(skipna=True).item()
                maxval = var.max(skipna=True).item()
                if 'minval' in var.attrs:
                    minval = np.minimum(minval, var.minval)
                    maxval = np.maximum(maxval, var.maxval)
                # Attributes shall have same datatype as in output file
                if 'dtype' in var.attrs:
                    att_type = var.attrs['dtype']
                else:
                    att_type = np.float32
                minval = att_type(minval)
                maxval = att_type(maxval)
                self.result[varname] = self.result[varname].assign_attrs({'minval': minval, 'maxval': maxval})

        if final is True:
            numtimes_before = self.result.sizes['time']
            self.result = self.result.sel({'time': slice(self.result.time[0], pd.to_datetime(self.time))})
            numtimes_after = self.result.sizes['time']
            if numtimes_after < numtimes_before:
                logger.debug(f'Truncating buffer from {numtimes_before} to {numtimes_after} times')

            # Final update some variable attributes
            status_dtype = self.ElementType.variables['status']['dtype']
            self.result['status'] = self.result.status.assign_attrs(
                {'valid_range': np.array((0, len(self.status_categories) - 1)).astype(status_dtype),
                 'flag_values': np.array(np.arange(len(self.status_categories))).astype(status_dtype),
                 'flag_meanings': " ".join(self.status_categories)
                })

            # Add items from metadata_dict
            for var in self.required_variables:
                keyword = 'reader_' + var
                if var not in self.env.priority_list:
                    fallback = self.get_config(f'environment:fallback:{var}')
                    if fallback is not None:
                        self.add_metadata(keyword, fallback)
                    else:
                        self.add_metadata(keyword, None)
                else:
                    readers = self.env.priority_list[var]
                    if readers[0].startswith(
                            'constant_reader') and var in self.env.readers[
                                readers[0]]._parameter_value_map:
                        self.add_metadata(
                            keyword, self.env.readers[
                                readers[0]]._parameter_value_map[var][0])
                    else:
                        self.add_metadata(keyword, self.env.priority_list[var])
            self.metadata_dict = {key:str(value) for key,value in self.metadata_dict.items()}

            final_metadata = {
                'time_coverage_end': str(self.time),
                'time_coverage_duration': pd.Timedelta(self.time-self.start_time).isoformat(),
                'time_coverage_resolution': pd.Timedelta(self.time_step).isoformat(),
                'performance': self.performance(),
                'geospatial_bounds_crs': 'EPSG:4326',
                'geospatial_bounds_vertical_crs': 'EPSG:5831',
                'geospatial_lat_min': self.result.lat.minval,
                'geospatial_lat_max': self.result.lat.maxval,
                'geospatial_lat_units': 'degrees_north',
                'geospatial_lat_resolution': 'point',
                'geospatial_lon_min': self.result.lon.minval,
                'geospatial_lon_max': self.result.lon.maxval,
                'geospatial_lon_units': 'degrees_east',
                'geospatial_lon_resolution': 'point',
                'runtime': str(self.timing['total time'])
                }
            if 'z' in self.result.data_vars:
                final_metadata['geospatial_vertical_min'] = self.result.z.minval
                final_metadata['geospatial_vertical_max'] = self.result.z.maxval
                final_metadata['geospatial_vertical_positive'] = 'up'

            self.result = self.result.assign_attrs(self.metadata_dict | final_metadata)

        if final is True or buffer_is_full:
            logger.debug('Writing to file')
            if self.outfile is not None:
                self.io_write_buffer()

        if final is False and buffer_is_full:
            logger.debug(f'Initialising new buffer')
            for var in self.result.data_vars:
                if 'time' in self.result[var].dims:
                    self.result[var][:] = np.nan
            newtime = self.result.coords['time'] + pd.Timedelta(self.time_step_output)*self.export_buffer_length
            self.result.coords['time'] = newtime
            logger.debug(f'Reset self.result, size {self.result.sizes}')

    def report_missing_variables(self):
        """Issue warning if some environment variables missing."""

        missing_variables = []
        for var in self.required_variables:
            if np.isnan(getattr(self.environment, var).min()):
                missing_variables.append(var)

        if len(missing_variables) > 0:
            logger.warning('Missing variables: ' + str(missing_variables))
            self.store_message('Missing variables: ' + str(missing_variables))

    def index_of_first_and_last(self):
        """Return the indices when elements were seeded and deactivated."""

        if self.index_of_first is None:
            if self.result is not None:
                index_of_first = self.result.lon.notnull().argmax(dim='time')
                index_of_last = len(self.result.time) - 1 - self.result.lon.isel(
                                                time=slice(None, None, -1)).notnull().argmax(dim='time')
            else:
                index_of_first = 0
                index_of_last = 0
            # Store for later retrieval
            self.index_of_first = index_of_first
            self.index_of_last = index_of_last
        else:  # Retrieved from cache
            index_of_first = self.index_of_first
            index_of_last = self.index_of_last

        return index_of_first, index_of_last

    def set_up_map(self,
                   corners=None,
                   buffer=.1,
                   delta_lat=None,
                   lscale=None,
                   fast=False,
                   hide_landmask=False,
                   xlocs = None,
                   ylocs = None,
                   **kwargs):
        """
        Generate Figure instance on which trajectories are plotted.

        :param hide_landmask: do not plot landmask (default False)
        :type hide_landmask: bool

        provide corners=[lonmin, lonmax, latmin, latmax] for specific map selection
        """
        logger.debug(f"Setting up map: {corners=}, {fast=}, {lscale=}")

        if self.result is None:  # To allow plotting of scheduled elements, if simulation not done
            lons = np.reshape(self.elements_scheduled.lon, (-1, 1))
            lats = np.reshape(self.elements_scheduled.lat, (-1, 1))
        else:
            lons = self.result.lon
            lats = self.result.lat

        if corners is not None:  # User provided map corners
            lonmin = corners[0]
            lonmax = corners[1]
            latmin = corners[2]
            latmax = corners[3]
        else:
            if 'compare_lonmin' in kwargs:  # checking min/max lon/lat of other simulations
                lonmin = np.minimum(kwargs['compare_lonmin'], np.nanmin(lons))
                lonmax = np.maximum(kwargs['compare_lonmax'], np.nanmax(lons))
                latmin = np.minimum(kwargs['compare_latmin'], np.nanmin(lats))
                latmax = np.maximum(kwargs['compare_latmax'], np.nanmax(lats))
            else:
                lonmin = np.nanmin(lons)
                lonmax = np.nanmax(lons)
                latmin = np.nanmin(lats)
                latmax = np.nanmax(lats)
            lonmin = lonmin - buffer * 2
            lonmax = lonmax + buffer * 2
            latmin = latmin - buffer
            latmax = latmax + buffer

        if fast is True:
            logger.warning(
                'Plotting fast. This will make your plots less accurate.')

            import matplotlib.style as mplstyle
            mplstyle.use(['fast'])

            # use a spherical earth
            axis = 57.29577951308232  # something to do with pi
            globe = ccrs.Globe(ellipse=None,
                               semimajor_axis=axis,
                               semiminor_axis=axis)
            self.crs_plot = ccrs.Mercator(globe=globe)
            self.crs_lonlat = ccrs.PlateCarree(globe=globe)

            if lscale is None:
                lscale = 'c'
        else:
            self.crs_plot = ccrs.Mercator()
            self.crs_lonlat = ccrs.PlateCarree()

            if lscale is None:
                lscale = 'auto'

        meanlat = (latmin + latmax) / 2
        aspect_ratio = float(latmax - latmin) / (float(lonmax - lonmin))
        aspect_ratio = aspect_ratio / np.cos(np.radians(meanlat))
        if 'figsize' in kwargs:
            figsize = kwargs['figsize']
        else:
            figsize = 11.  # inches

        if aspect_ratio > 1:
            fig = plt.figure(figsize=(figsize / aspect_ratio, figsize))
        else:
            fig = plt.figure(figsize=(figsize, figsize * aspect_ratio))

        ax = fig.add_subplot(111, projection=self.crs_plot)
        if lonmin == -180 and lonmax == 180:
            lonmax -= .1  # To avoid problem with Cartopy
        ax.set_extent([lonmin, lonmax, latmin, latmax], crs=self.crs_lonlat)

        gl = ax.gridlines(self.crs_lonlat, draw_labels=True, xlocs=xlocs, ylocs=ylocs)
        gl.top_labels = None

        if 'ocean_color' in kwargs:
            ax.patch.set_facecolor(kwargs['ocean_color'])
            ocean_color = kwargs['ocean_color']
        else:
            ocean_color = 'white'
        if 'land_color' in kwargs:
            land_color = kwargs['land_color']
        else:
            if fast is True:
                land_color = 'gray'
            else:
                land_color = cfeature.COLORS['land']

        if 'text' in kwargs:
            if not isinstance(kwargs['text'], list):
                text = list(kwargs['text'])
            else:
                text = kwargs['text']
            for te in text:
                plt.text(transform=self.crs_lonlat, **te)

        if 'box' in kwargs:
            if not isinstance(kwargs['box'], list):
                box = list(kwargs['box'])
            else:
                box = kwargs['box']
            for bx in box:
                lonmn = bx['lon'][0]
                lonmx = bx['lon'][1]
                latmn = bx['lat'][0]
                latmx = bx['lat'][1]
                del bx['lon']
                del bx['lat']
                if 'text' in bx:
                    plt.text(x=lonmn,
                             y=latmx,
                             s=bx['text'],
                             transform=self.crs_lonlat)
                    del bx['text']
                patch = matplotlib.patches.Rectangle(
                    xy=[lonmn, latmn],
                    width=lonmx - lonmn,
                    height=latmx - latmn,
                    transform=self.crs_lonlat,
                    zorder=10,
                    **bx)
                ax.add_patch(patch)

        if not hide_landmask:
            if 'land_binary_mask' in self.env.priority_list and self.env.priority_list[
                    'land_binary_mask'][0] == 'shape':
                logger.debug('Using custom shapes for plotting land..')
                ax.add_geometries(self.env.readers['shape'].polys,
                                  self.crs_lonlat,
                                  facecolor=land_color,
                                  edgecolor='black')
            else:
                reader_global_landmask.plot_land(ax, lonmin, latmin, lonmax,
                                                 latmax, fast, ocean_color,
                                                 land_color, lscale,
                                                 crs_plot=self.crs_plot,
                                                 crs_lonlat=self.crs_lonlat)

        fig.canvas.draw()
        fig.set_layout_engine('tight')

        index_of_first, index_of_last = self.index_of_first_and_last()

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except:
            pass

        try:  # Maximise figure window size
            mng.resize(*mng.window.maxsize())
        except:
            pass

        # TODO: avoid transposing lon, lat, and avoid returning lon, lat in the first place
        return fig, ax, self.crs_plot, lons.T, lats.T, index_of_first, index_of_last

    # Obsolete, to be removed
    #def get_lonlats(self):
    #    if self.result is not None:
    #        lons = self.result.lon
    #        lats = self.result.lat
    #    else:
    #        #if len(self.result.time) > 0:
    #        #    lons = np.ma.array(np.reshape(self.elements.lon, (1, -1))).T
    #        #    lats = np.ma.array(np.reshape(self.elements.lat, (1, -1))).T
    #        #else:
    #        lons = np.ma.array(
    #            np.reshape(self.elements_scheduled.lon, (1, -1))).T
    #        lats = np.ma.array(
    #            np.reshape(self.elements_scheduled.lat, (1, -1))).T
    #    return lons, lats

    def animation(self,
                  buffer=.2,
                  corners=None,
                  filename=None,
                  compare=None,
                  compare_marker='o',
                  background=None,
                  alpha=1,
                  bgalpha=.5,
                  vmin=None,
                  vmax=None,
                  drifter=None,
                  shapefiles=None,
                  skip=None,
                  scale=None,
                  color=False,
                  clabel=None,
                  colorbar=True,
                  cmap=None,
                  density=False,
                  show_elements=True,
                  show_trajectories=False,
                  linewidth=1,
                  trajectory_alpha=.1,
                  hide_landmask=False,
                  density_pixelsize_m=1000,
                  unitfactor=1,
                  lcs=None,
                  surface_only=False,
                  markersize=20,
                  markersize_scaling=None,
                  origin_marker=None,
                  legend=None,
                  legend_loc='best',
                  title='auto',
                  fps=8,
                  lscale=None,
                  fast=False,
                  blit=False,
                  frames=None,
                  xlocs = None,
                  ylocs = None,
                  **kwargs):
        """Animate last run."""

        filename = str(filename) if filename is not None else None

        if self.result is not None and self.num_elements_total(
        ) == 0 and not hasattr(self, 'ds'):
            raise ValueError('Please run simulation before animating')

        if compare is not None:
            compare_list, compare_args = self._get_comparison_xy_for_plots(
                compare)
            kwargs.update(compare_args)

        start_time = datetime.now()
        if cmap is None:
            cmap = 'jet'
        if isinstance(cmap, str):
            cmap = matplotlib.colormaps[cmap]

        if color is False and background is None and lcs is None and density is False:
            colorbar = False

        markercolor = self.plot_comparison_colors[0]

        if isinstance(density, str):
            # Density field is weighted by this variable
            # TODO: not yet implemented!
            density_weight = density
            density = True
        else:
            if density is True:
                density_weight = None
            elif density is not False:
                density_weight = density
                density = True
        if density is True:  # Get density arrays
            if hasattr(self, 'ds'):  # opened with Xarray
                if origin_marker is None:
                    origin_marker = 0
                    per_origin_marker = False
                else:
                    per_origin_marker = True
                H, H_om, lon_array, lat_array = self.get_density_xarray(
                    pixelsize_m=density_pixelsize_m, weights=density_weight)
                if per_origin_marker is True:
                    H = H_om[:, :, :, origin_marker]
            else:
                if origin_marker is not None:
                    raise ValueError(
                        'Separation by origin_marker is only active when imported from file with '
                        'open_xarray: https://opendrift.github.io/gallery/example_huge_output.html'
                    )
                H, H_submerged, H_stranded, lon_array, lat_array = \
                    self.get_density_array(pixelsize_m=density_pixelsize_m,
                                           weight=density_weight)
                H = H + H_submerged + H_stranded

        # Find map coordinates and plot points with empty data
        fig, ax, crs, x, y, index_of_first, index_of_last = \
            self.set_up_map(buffer=buffer, corners=corners, lscale=lscale,
                            fast=fast, hide_landmask=hide_landmask, xlocs = xlocs, ylocs = ylocs, **kwargs)

        def plot_timestep(i):
            """Sub function needed for matplotlib animation."""

            time_string = np.datetime_as_string(self.result.time[i], unit='s')
            ret = [points, points_deactivated
                   ]  # list of elements to return for blitting
            if title == 'auto':
                ax.set_title('%s\n%s UTC' % (self._figure_title(), time_string))
            else:
                ax.set_title('%s\n%s UTC' % (title, time_string))
            if background is not None:
                ret.append(bg)
                if isinstance(background, xr.DataArray):
                    scalar = background[i, :, :].values
                else:
                    map_x, map_y, scalar, u_component, v_component = \
                        self.get_map_background(ax, background, self.crs_plot,
                                                time=self.result.time[i])
                # https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data
                bg.set_array(scalar.ravel())
                if type(background) is list:
                    ret.append(bg_quiv)
                    bg_quiv.set_UVC(u_component[::skip, ::skip],
                                    v_component[::skip, ::skip])

            if lcs is not None:
                ax.pcolormesh(lcs['lon'],
                              lcs['lat'],
                              lcs['ALCS'][i, :, :],
                              alpha=bgalpha,
                              vmin=vmin,
                              vmax=vmax,
                              cmap=cmap,
                              transform=self.crs_lonlat)

            if density is True:
                # Update density plot
                pm.set_array(H[i, :, :].ravel())
                ret.append(pm)

            # Move points
            if show_elements is True:
                points.set_offsets(np.c_[x[i, range(x.shape[1])],
                                         y[i, range(x.shape[1])]])
                points_deactivated.set_offsets(
                    np.c_[x_deactive[index_of_last_deactivated < i],
                          y_deactive[index_of_last_deactivated < i]])

                if isinstance(markersize, str):
                    points.set_sizes(markersize_scaling * np.abs(self.result[markersize][:, i]))

                if color is not False:  # Update colors
                    points.set_array(colorarray[:, i])
                    if compare is not None:
                        for cd in compare_list:
                            cd['points_other'].set_array(colorarray[:, i])
                    if (isinstance(color, str) or hasattr(color, '__len__')) and len(index_of_last_deactivated)>0:
                        points_deactivated.set_array(colorarray_deactivated[
                            index_of_last_deactivated < i])

            if drifter is not None:
                for drnum, dr in enumerate(drifter):
                    drifter_pos[drnum].set_offsets(np.c_[dr['x'][i],
                                                         dr['y'][i]])
                    drifter_line[drnum].set_data(dr['x'][0:i+1], dr['y'][0:i+1])
                    ret.append(drifter_line[drnum])
                    ret.append(drifter_pos[drnum])

            if shapefiles is not None:
                import geopandas as gpd
                for sf in shapefiles:
                    shdf = gpd.read_file(sf)
                    shdf = shdf.to_crs("EPSG:4326")
                    ax.add_geometries(shdf.geometry, self.crs_lonlat, edgecolor='g', linewidth=2, facecolor='none')

            if show_elements is True:
                if compare is not None:
                    for cd in compare_list:
                        cd['points_other'].set_offsets(
                            np.c_[cd['x_other'][range(cd['x_other'].shape[0]),
                                                i],
                                  cd['y_other'][range(cd['x_other'].shape[0]),
                                                i]])
                        cd['points_other_deactivated'].set_offsets(np.c_[
                            cd['x_other_deactive'][
                                cd['index_of_last_deactivated_other'] < i],
                            cd['y_other_deactive'][
                                cd['index_of_last_deactivated_other'] < i]])
                        ret.append(cd['points_other'])
                        ret.append(cd['points_other_deactivated'])

            return ret

        if surface_only is True:
            z = self.result.z.T
            x = x.where(z==0)
            y = y.where(z==0)

        if show_trajectories is True:
            ax.plot(x, y, color='gray', alpha=trajectory_alpha, transform=self.crs_lonlat, linewidth=linewidth)

        if color is not False and show_elements is True:
            if isinstance(color, str):
                colorarray = self.result[color]
                colorarray = colorarray * unitfactor
                colorarray_deactivated = \
                    self.result[color].T[
                        index_of_last[self.elements_deactivated.ID],
                                      self.elements_deactivated.ID]
            elif hasattr(color,
                         '__len__'):  # E.g. array/list of ensemble numbers
                colorarray_deactivated = color[self.elements_deactivated.ID]
                colorarray = np.tile(color, (len(self.result.time), 1)).T
            else:
                colorarray = color
            if vmin is None:
                vmin = colorarray.min()
                vmax = colorarray.max()

        if background is not None:
            if isinstance(background, xr.DataArray):
                map_x = background.coords['lon_bin']
                map_y = background.coords['lat_bin']
                scalar = background[0, :, :]
                map_y, map_x = np.meshgrid(map_y, map_x)
            else:
                map_x, map_y, scalar, u_component, v_component = \
                    self.get_map_background(ax, background, self.crs_plot,
                                            time=self.result.time[0])
            bg = ax.pcolormesh(map_x,
                               map_y,
                               scalar,
                               alpha=bgalpha,
                               zorder=1,
                               antialiased=True,
                               linewidth=0.0,
                               rasterized=True,
                               vmin=vmin,
                               vmax=vmax,
                               cmap=cmap,
                               transform=self.crs_lonlat)
            if type(background) is list:
                bg_quiv = ax.quiver(map_x[::skip, ::skip],
                                    map_y[::skip, ::skip],
                                    u_component[::skip, ::skip],
                                    v_component[::skip, ::skip],
                                    scale=scale,
                                    zorder=1,
                                    transform=self.crs_lonlat)

        if lcs is not None:
            if vmin is None:
                vmin = lcs['ALCS'].min()
                vmax = lcs['ALCS'].max()
            lcsh = ax.pcolormesh(lcs['lon'],
                                 lcs['lat'],
                                 lcs['ALCS'][0, :, :],
                                 vmin=vmin,
                                 vmax=vmax,
                                 cmap=cmap,
                                 transform=self.crs_lonlat)

        if show_elements is True:
            index_of_last_deactivated = \
                index_of_last[self.elements_deactivated.ID]
        if legend is None:
            legend = ['']

        if color is False:
            cargs = {'c': None, 'color': markercolor, 'cmap': None}
        else:
            cargs = {'c': [], 'color': None, 'cmap': cmap}

        if isinstance(markersize, str):
            if markersize_scaling is None:
                markersize_scaling = 20
            markersize_scaling = markersize_scaling / np.abs(self.result[markersize]).max()

        if isinstance(markersize, str):
            points = ax.scatter([], [],
                                **cargs,
                                zorder=10,
                                edgecolor=[],
                                alpha=alpha,
                                vmin=vmin,
                                vmax=vmax,
                                label=legend[0],
                                transform=self.crs_lonlat)
        else:
            points = ax.scatter([], [],
                                **cargs,
                                zorder=10,
                                edgecolor=[],
                                alpha=alpha,
                                s=markersize,
                                vmin=vmin,
                                vmax=vmax,
                                label=legend[0],
                                transform=self.crs_lonlat)

        if (compare is None) and (legend != ['']):
            markers = []
            for legend_index in np.arange(len(legend)):
                if legend[legend_index] != '':
                    markers.append(
                        matplotlib.lines.Line2D(
                            [0], [0],
                            marker='o',
                            color='w',
                            linewidth=0,
                            markeredgewidth=0,
                            markerfacecolor=cmap(legend_index /
                                                 (len(legend) - 1)),
                            markersize=10,
                            label=legend[legend_index]))
            legend = list(filter(None, legend))
            ax.legend(markers, legend, loc=legend_loc)

        # Plot deactivated elements, with transparency
        if isinstance(markersize, str):
            points_deactivated = ax.scatter([], [],
                                            **cargs,
                                            zorder=9,
                                            vmin=vmin,
                                            vmax=vmax,
                                            s=markersize_scaling,
                                            edgecolor=[],
                                            alpha=0,
                                            transform=self.crs_lonlat)
        else:
            points_deactivated = ax.scatter([], [],
                                            **cargs,
                                            zorder=9,
                                            vmin=vmin,
                                            vmax=vmax,
                                            s=markersize,
                                            edgecolor=[],
                                            alpha=.3,
                                            transform=self.crs_lonlat)

        x_deactive, y_deactive = (self.elements_deactivated.lon,
                                  self.elements_deactivated.lat)

        if compare is not None:
            for cn, cd in enumerate(compare_list):
                if legend != ['']:
                    legstr = legend[cn + 1]
                else:
                    legstr = None
                if color is False:
                    cargs['color'] = self.plot_comparison_colors[cn + 1]
                #else:
                #    c = []
                cd['points_other'] = \
                    ax.scatter([], [], **cargs, marker=compare_marker,
                               s=markersize, label=legstr, zorder=10, transform = self.crs_lonlat)
                # Plot deactivated elements, with transparency
                cd['points_other_deactivated'] = \
                    ax.scatter([], [], **cargs, alpha=.3, zorder=9, marker=compare_marker,
                               s=markersize, transform = self.crs_lonlat)

            if legend != ['', '']:
                plt.legend(markerscale=2, loc=legend_loc)

        if density is True:
            cmap.set_under('w')
            H = np.ma.masked_where(H == 0, H)
            lat_array, lon_array = np.meshgrid(lat_array, lon_array)
            if vmax is None:
                vmax = H.max()
            pm = ax.pcolormesh(lon_array,
                               lat_array,
                               H[0, :, :],
                               vmin=0.1,
                               vmax=vmax,
                               cmap=cmap,
                               transform=self.crs_lonlat)

        if drifter is not None:
            if not isinstance(drifter, list):
                drifter = [drifter]
            drifter_pos = [None] * len(drifter)
            drifter_line = [None] * len(drifter)
            for drnum, dr in enumerate(drifter):
                # Interpolate drifter time series onto simulation times
                sts = (self.result.time - self.result.time[0]) / np.timedelta64(1, 's')
                dr_times = np.array(dr['time'], dtype=self.result.time.dtype)
                dts = (dr_times-dr_times[0]) / np.timedelta64(1, 's')
                dr['x'] = np.interp(sts, dts, dr['lon'])
                dr['y'] = np.interp(sts, dts, dr['lat'])
                dr['x'][sts < dts[0]] = np.nan
                dr['x'][sts > dts[-1]] = np.nan
                dr['y'][sts < dts[0]] = np.nan
                dr['y'][sts > dts[-1]] = np.nan
                dlabel = dr['label'] if 'label' in dr else 'Drifter'
                dcolor = dr['color'] if 'color' in dr else 'r'
                dlinewidth = dr['linewidth'] if 'linewidth' in dr else 2
                dzorder = dr['zorder'] if 'zorder' in dr else 10
                dmarkersize = dr['markersize'] if 'markersize' in dr else 20
                drifter_pos[drnum] = ax.scatter([], [],
                                                c=dcolor,
                                                zorder=dzorder + 1,
                                                s=dmarkersize,
                                                label=dlabel,
                                                transform=self.crs_lonlat)
                drifter_line[drnum] = ax.plot([], [],
                                              color=dcolor,
                                              linewidth=dlinewidth,
                                              zorder=dzorder,
                                              transform=self.crs_lonlat)[0]
            plt.legend()

        fig.canvas.draw()
        fig.set_layout_engine('tight')
        if colorbar is True:
            if color is not False:
                if isinstance(color, str) or clabel is not None:
                    if clabel is None:
                        clabel = color
                item = points
            elif density is not False:
                item = pm
                if clabel is None:
                    clabel = 'density'
            elif lcs is not None:
                item = lcsh
                if clabel is None:
                    clabel = 'LCS'
            elif background is not None:
                item = bg
                if clabel is None:
                    if isinstance(background, xr.DataArray):
                        clabel = background.name
                    else:
                        clabel = background

            cb = fig.colorbar(item,
                              orientation='horizontal',
                              pad=.05,
                              aspect=30,
                              shrink=.8,
                              drawedges=False)
            cb.set_label(clabel)

        frames = x.shape[0] if frames is None else frames

        if compare is not None:
            frames = min(x.shape[0], cd['x_other'].shape[1])

        # blit is now provided to animation()
        #blit = sys.platform != 'darwin'  # blitting does not work on mac os

        self.__save_or_plot_animation__(plt.gcf(),
                                        plot_timestep,
                                        filename,
                                        frames,
                                        fps,
                                        interval=50,
                                        blit=blit)

        logger.info('Time to make animation: %s' %
                    (datetime.now() - start_time))

    def __save_or_plot_animation__(self, figure, plot_timestep, filename,
                                   frames, fps, interval, blit):

        if filename is not None or 'sphinx_gallery' in sys.modules:
            logger.debug("Saving animation..")
            self.__save_animation__(figure,
                                    plot_timestep,
                                    filename,
                                    frames=frames,
                                    fps=fps,
                                    blit=blit,
                                    interval=interval)

        else:
            logger.debug("Showing animation..")
            anim = animation.FuncAnimation(figure,
                                           plot_timestep,
                                           blit=blit,
                                           frames=frames,
                                           interval=interval)
            try:
                plt.show()
            except AttributeError as e:
                logger.exception(e)
                pass

    @require_mode(mode=Mode.Result)
    def animation_profile(self,
                          filename=None,
                          compare=None,
                          markersize=20,
                          markersize_scaling=None,
                          alpha=1,
                          fps=20,
                          color=None,
                          cmap=None,
                          vmin=None,
                          vmax=None,
                          legend=None,
                          legend_loc=None):
        """Animate vertical profile of the last run."""

        start_time = datetime.now()

        def plot_timestep(i):
            """Sub function needed for matplotlib animation."""
            time_string = np.datetime_as_string(self.result.time[i], unit='s')
            ax.set_title('%s UTC' % time_string)
            points.set_offsets(np.c_[x[range(x.shape[0]), i], z[range(x.shape[0]), i]])
            if color is not None and compare is None:
                points.set_array(colorarray[:, i])

            points_deactivated.set_offsets(np.c_[
                x_deactive[index_of_last_deactivated < i],
                z_deactive[index_of_last_deactivated < i]])

            if isinstance(markersize, str):
                points.set_sizes(np.abs(markersize_scaling * self.result[markersize][:, i]))

            if compare is not None:
                points_other.set_offsets(np.c_[x_other[range(x_other.shape[0]), i],
                                               z_other[range(x_other.shape[0]), i]])
                points_other_deactivated.set_offsets(np.c_[
                     x_other_deactive[index_of_last_deactivated_other < i],
                     z_other_deactive[index_of_last_deactivated_other < i]])
                return points, points_deactivated, points_other,
            else:
                return points, points_deactivated,

        if cmap is None:
            cmap = 'jet'
        if isinstance(cmap, str):
            cmap = matplotlib.colormaps[cmap]
        if color is not False:
            if isinstance(color, str):
                colorarray = self.result[color]
                if vmin is None:
                    vmin = colorarray.min()
                    vmax = colorarray.max()

        markercolor = self.plot_comparison_colors[0]

        if color is None:
            cargs = {'c': None, 'color': markercolor, 'cmap': None}
        else:
            cargs = {'c': [], 'color': None, 'cmap': cmap}

        # Set up plot
        index_of_first, index_of_last = \
            self.index_of_first_and_last()
        z = self.result.z
        x = self.result.lon

        #seafloor_depth = \
        #    -self.result.sea_floor_depth_below_sea_level
        fig = plt.figure(figsize=(10, 6.))  # Suitable aspect ratio
        ax = fig.gca()
        plt.xlabel('Longitude [degrees]')
        plt.ylabel('Depth [m]')
        index_of_last_deactivated = \
            index_of_last[self.elements_deactivated.ID]
        if legend is None:
            if compare is None:
                legs = ['', '']
        else:
            legs = legend
        if isinstance(markersize, str):
            ms = None
        else:
            ms = markersize
        if isinstance(markersize, str):
            if markersize_scaling is None:
                markersize_scaling = 20
            markersize_scaling = markersize_scaling / np.abs(self.result[markersize]).max()

        points = ax.scatter([], [],
                            **cargs,
                            zorder=10,
                            edgecolor=[],
                            alpha=alpha,
                            s=ms,
                            label=legs[0],
                            vmin=vmin,
                            vmax=vmax)

        markers = []
        if compare is None and legend is not None:  # Empty points to get legend
            for legend_index in np.arange(len(legend)):
                if legend[legend_index] != '':
                    markers.append(
                        matplotlib.lines.Line2D(
                            [0], [0], marker='o', linewidth=0, markeredgewidth=0,
                        markerfacecolor=cmap(legend_index /
                                             (len(legend) - 1)),
                        markersize=10,
                        label=legend[legend_index]))

            legend = list(filter(None, legend))
            leg = ax.legend(markers, legend, loc=legend_loc)
            leg.set_zorder(20)

        # Plot deactivated elements, with transparency
        points_deactivated = ax.scatter([], [], **cargs, zorder=10, edgecolor=[],
                                        s=ms, vmin=vmin, vmax=vmax)
        x_deactive = self.elements_deactivated.lon
        z_deactive = self.elements_deactivated.z

        if compare is not None:
            if type(compare) is str:
                # Other is given as filename
                other = self.__class__()
                other.io_import_file(compare)
            else:
                # Other is given as an OpenDrift object
                other = compare
            z_other = other.result.z
            x_other = other.result.lon
            points_other = ax.scatter([], [],
                            color='r',
                            zorder=10,
                            alpha=alpha,
                            edgecolor=[],
                            s=markersize,
                            label=legs[1],
                            vmin=vmin,
                            vmax=vmax)

            # Plot deactivated elements, with transparency
            points_other_deactivated = ax.scatter([], [], color='r', s=markersize, alpha=.3)
            x_other_deactive = other.elements_deactivated.lon
            z_other_deactive = other.elements_deactivated.z
            index_of_first_other, index_of_last_other = other.index_of_first_and_last()
            index_of_last_deactivated_other = \
                index_of_last_other[other.elements_deactivated.ID]
            xmax = np.maximum(x.max(), x_other.max())
            xmin = np.minimum(x.min(), x_other.min())
            zmax = np.maximum(z.max(), z_other.max())
            zmin = np.minimum(z.min(), z_other.min())
        else:
            xmin = x.min()
            xmax = x.max()
            zmin = z.min()
            zmax = z.max()

        # Set figure limits
        sky = (zmax - zmin) * .1  # Sky height is 10% of water depth
        plt.xlim([xmin, xmax])
        plt.ylim([zmin, sky])
        ax.add_patch(
            plt.Rectangle((xmin, 0), xmax - xmin, sky, color='lightsteelblue'))
        ax.add_patch(
            plt.Rectangle((xmin, zmin),
                          xmax - xmin,
                          -zmin,
                          color='cornflowerblue'))

        if legend is not None and compare is not None:
            plt.legend(loc=4)

        self.__save_or_plot_animation__(plt.gcf(),
                                        plot_timestep,
                                        filename,
                                        x.shape[1],
                                        fps,
                                        interval=150,
                                        blit=False)

        logger.info('Time to make animation: %s' %
                    (datetime.now() - start_time))

    def _get_comparison_xy_for_plots(self, compare):
        if not type(compare) is list:
            compare = [compare]
        compare_list = [{}] * len(compare)
        lonmin = 1000
        lonmax = -1000
        latmin = 1000
        latmax = -1000
        for cn, comp in enumerate(compare):
            compare_list[cn] = {}
            cd = compare_list[cn]  # pointer to dict with data
            if type(comp) is str:
                # Other is given as filename
                other = self.__class__()
                other.io_import_file(comp)
            else:
                # Other is given as an OpenDrift object
                other = comp

            lonmin = np.minimum(lonmin, other.result.lon.min())
            lonmax = np.maximum(lonmax, other.result.lon.max())
            latmin = np.minimum(latmin, other.result.lat.min())
            latmax = np.maximum(latmax, other.result.lat.max())

            # Find map coordinates of comparison simulations
            cd['x_other'] = other.result.lon.copy()
            cd['y_other'] = other.result.lat.copy()
            cd['x_other_deactive'], cd['y_other_deactive'] = \
                (other.elements_deactivated.lon.copy(),
                    other.elements_deactivated.lat.copy())
            cd['index_of_last_other'] = other.index_of_first_and_last()[1]
            cd['index_of_last_deactivated_other'] = \
                cd['index_of_last_other'][other.elements_deactivated.ID]

        compare_args = {
            'compare_lonmin': lonmin,
            'compare_lonmax': lonmax,
            'compare_latmin': latmin,
            'compare_latmax': latmax
        }

        return compare_list, compare_args

    def plot(self,
             background=None,
             buffer=.2,
             corners=None,
             linecolor=None,
             filename=None,
             show=True,
             vmin=None,
             vmax=None,
             compare=None,
             cmap='jet',
             lvmin=None,
             lvmax=None,
             skip=None,
             scale=None,
             show_scalar=True,
             contourlines=False,
             drifter=None,
             colorbar=True,
             linewidth=1,
             lcs=None,
             show_elements=True,
             show_trajectories=True,
             show_initial=True,
             density_pixelsize_m=1000,
             lalpha=None,
             bgalpha=1,
             clabel=None,
             cpad=.05,
             caspect=30,
             cshrink=.8,
             surface_color=None,
             submerged_color=None,
             markersize=20,
             title='auto',
             legend=True,
             legend_loc='best',
             lscale=None,
             fast=False,
             hide_landmask=False,
             xlocs = None,
             ylocs = None,
             **kwargs):
        """Basic built-in plotting function intended for developing/debugging.

        Plots trajectories of all particles.
        Positions marked with colored stars:
        - green: all start positions
        - red: deactivated particles
        - blue: particles still active at end of simulation

        Requires availability of Cartopy.

        Arguments:
            background: string, name of variable (standard_name) which will
                be plotted as background of trajectories, provided that it
                can be read with one of the available readers.

            buffer: float; spatial buffer of plot in degrees of
                longitude/latitude around particle collection.

            background: name of variable to be plotted as background field.
            Use two element list for vector fields, e.g. ['x_wind', 'y_wind']

            vmin, vmax: minimum and maximum values for colors of background.

            linecolor: name of variable to be used for coloring trajectories, or matplotlib color string.

            lvmin, lvmax: minimum and maximum values for colors of trajectories.

            lscale (string): resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto'). default is 'auto'.

            fast (bool): use some optimizations to speed up plotting at the cost of accuracy

            :param hide_landmask: do not plot landmask (default False).
            :type hide_landmask: bool
        """

        mappable = None

        if self.result is not None and self.num_elements_total(
        ) == 0 and not hasattr(self, 'ds'):
            raise ValueError('Please run simulation before animating')

        start_time = datetime.now()

        if compare is not None:
            # Extend map coverage to cover comparison simulations
            cd, compare_args = self._get_comparison_xy_for_plots(compare)
            kwargs.update(compare_args)

        if drifter is not None:
            # Extend map coverage to cover provided trajectory
            # TODO: drifter should be list of dictionaries
            ttime = np.array(drifter['time'], dtype=self.result.time.dtype)
            i = np.where((ttime >= self.result.time[0].values) & (ttime <= self.result.time[-1].values))[0]
            drifter['lon'] = np.atleast_1d(drifter['lon'])
            drifter['lat'] = np.atleast_1d(drifter['lat'])
            tlonmin = drifter['lon'][i].min()
            tlonmax = drifter['lon'][i].max()
            tlatmin = drifter['lat'][i].min()
            tlatmax = drifter['lat'][i].max()
            if 'compare_lonmin' not in kwargs:
                kwargs['compare_lonmin'] = tlonmin
                kwargs['compare_lonmax'] = tlonmax
                kwargs['compare_latmin'] = tlatmin
                kwargs['compare_latmax'] = tlatmax
            else:
                kwargs['compare_lonmin'] = np.minimum(kwargs['compare_lonmin'],
                                                      tlonmin)
                kwargs['compare_lonmax'] = np.maximum(kwargs['compare_lonmax'],
                                                      tlonmax)
                kwargs['compare_latmin'] = np.minimum(kwargs['compare_latmin'],
                                                      tlatmin)
                kwargs['compare_latmax'] = np.maximum(kwargs['compare_latmax'],
                                                      tlatmax)

        fig, ax, crs, x, y, index_of_first, index_of_last = \
            self.set_up_map(buffer=buffer, corners=corners, lscale=lscale, fast=fast, hide_landmask=hide_landmask, xlocs = xlocs, ylocs = ylocs, **kwargs)

        markercolor = self.plot_comparison_colors[0]

        # The more elements, the more transparent we make the lines
        if lalpha is None:
            min_alpha = 0.1
            max_elements = 5000.0
            alpha = min_alpha**(2 * (self.num_elements_total() - 1) /
                                (max_elements - 1))
            alpha = np.max((min_alpha, alpha))
        else:
            alpha = lalpha  #  provided transparency of trajectories
        if legend is False:
            legend = None
        if self.result is not None and linewidth != 0 and show_trajectories is True:
            # Plot trajectories
            from matplotlib.colors import is_color_like
            if linecolor is None or is_color_like(linecolor) is True:
                if is_color_like(linecolor):
                    linecolor = linecolor
                else:
                    linecolor = 'gray'
                if compare is not None and legend is not None:
                    if legend is True:
                        if hasattr(compare, 'len'):
                            numleg = len(compare)
                        else:
                            numleg = 2
                        legend = [
                            'Simulation %d' % (i + 1) for i in range(numleg)
                        ]
                    ax.plot(x[:, 0],
                            y[:, 0],
                            color=linecolor,
                            alpha=alpha,
                            label=legend[0],
                            linewidth=linewidth,
                            transform=self.crs_lonlat)
                    ax.plot(x,
                            y,
                            color=linecolor,
                            alpha=alpha,
                            label='_nolegend_',
                            linewidth=linewidth,
                            transform=self.crs_lonlat)
                else:
                    ax.plot(x,
                            y,
                            color=linecolor,
                            alpha=alpha,
                            linewidth=linewidth,
                            transform=self.crs_lonlat)
            else:
                #colorbar = True
                # Color lines according to given parameter
                try:
                    if isinstance(linecolor, str):
                        param = self.result[linecolor]
                    elif hasattr(linecolor, '__len__'):
                        param = np.tile(linecolor, (len(self.result.time), 1)).T
                    else:
                        param = linecolor
                except:
                    raise ValueError(
                        'Available parameters to be used for linecolors: ' +
                        str(self.result.data_vars))
                from matplotlib.collections import LineCollection
                for i in range(x.shape[1]):
                    vind = np.arange(index_of_first[i], index_of_last[i] + 1)
                    points = np.array([x[vind, i].T,
                                       y[vind, i].T]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]],
                                              axis=1)
                    if lvmin is None:
                        lvmin = param.min()
                        lvmax = param.max()
                    lc = LineCollection(
                        segments,
                        #cmap=plt.colormaps['Spectral'],
                        cmap=cmap,
                        norm=plt.Normalize(lvmin, lvmax),
                        transform=self.crs_lonlat)
                    #lc.set_linewidth(3)
                    lc.set_array(param.T[vind, i])
                    mappable = ax.add_collection(lc)
                #axcb = fig.colorbar(lc, ax = ax, orientation = 'horizontal')
                #try:  # Add unit to colorbar if available
                #    colorbarstring = linecolor + '  [%s]' % \
                #        (self.history_metadata[linecolor]['units'])
                #except:
                #    colorbarstring = linecolor
                ##axcb.set_label(colorbarstring)
                #axcb.set_label(colorbarstring, size=14)
                #axcb.ax.tick_params(labelsize=14)

        if compare is None:
            label_initial = 'initial (%i)' % x.shape[1]
            label_active = 'active (%i)' % (x.shape[1] -
                                            self.num_elements_deactivated())
            color_initial = self.status_colors['initial']
            color_active = self.status_colors['active']
        else:
            label_initial = None
            label_active = None
            color_initial = 'gray'
            color_active = 'gray'
        if show_elements is True:
            if show_initial is True:
                ax.scatter(x[index_of_first, range(x.shape[1])],
                           y[index_of_first, range(x.shape[1])],
                           s=markersize,
                           zorder=10,
                           edgecolor=markercolor,
                           linewidths=.2,
                           color=color_initial,
                           label=label_initial,
                           transform=self.crs_lonlat)
            if surface_color is not None:
                color_active = surface_color
                label_active = 'surface'
            ax.scatter(x[index_of_last, range(x.shape[1])],
                       y[index_of_last, range(x.shape[1])],
                       s=markersize,
                       zorder=3,
                       edgecolor=markercolor,
                       linewidths=.2,
                       color=color_active,
                       label=label_active,
                       transform=self.crs_lonlat)
            #if submerged_color is not None:
            #    map.scatter(x[range(x.shape[0]), index_of_last],
            #                y[range(x.shape[0]), index_of_last], s=markersize,
            #                zorder=3, edgecolor=markercolor, linewidths=.2,
            #                c=submerged_color, label='submerged')

            x_deactivated, y_deactivated = (self.elements_deactivated.lon,
                                            self.elements_deactivated.lat)
            # Plot deactivated elements, labeled by deactivation reason
            for statusnum, status in enumerate(self.status_categories):
                if status == 'active':
                    continue  # plotted above
                if status not in self.status_colors:
                    # If no color specified, pick an unused one
                    for color in [
                            'red', 'blue', 'green', 'black', 'gray', 'cyan',
                            'DarkSeaGreen', 'brown'
                    ]:
                        if color not in self.status_colors.values():
                            self.status_colors[status] = color
                            break
                indices = np.where(
                    self.elements_deactivated.status == statusnum)
                if len(indices[0]) > 0:
                    if (status == 'seeded_on_land'
                            or status == 'seeded_at_nodata_position'):
                        zorder = 11
                    else:
                        zorder = 3
                    if compare is not None:
                        legstr = None
                    else:
                        legstr = '%s (%i)' % (status, len(indices[0]))
                    if compare is None:
                        color_status = self.status_colors[status]
                    else:
                        color_status = 'gray'
                    ax.scatter(x_deactivated[indices],
                               y_deactivated[indices],
                               s=markersize,
                               zorder=zorder,
                               edgecolor=markercolor,
                               linewidths=.1,
                               color=color_status,
                               label=legstr,
                               transform=self.crs_lonlat)

        if compare is not None:
            for i, c in enumerate(cd):
                if legend != None:
                    legstr = legend[i + 1]
                else:
                    legstr = None
                ax.plot(c['x_other'].T[:, 0],
                        c['y_other'].T[:, 0],
                        color=self.plot_comparison_colors[i + 1],
                        linestyle='-',
                        label=legstr,
                        linewidth=linewidth,
                        transform=self.crs_lonlat)
                ax.plot(c['x_other'].T,
                        c['y_other'].T,
                        color=self.plot_comparison_colors[i + 1],
                        linestyle='-',
                        label='_nolegend_',
                        linewidth=linewidth,
                        transform=self.crs_lonlat)
                ax.scatter(c['x_other'][range(c['x_other'].shape[0]),
                                        c['index_of_last_other']],
                           c['y_other'][range(c['y_other'].shape[0]),
                                        c['index_of_last_other']],
                           s=markersize,
                           zorder=3,
                           edgecolor=markercolor,
                           linewidths=.2,
                           color=self.plot_comparison_colors[i + 1],
                           transform=self.crs_lonlat)

        if background is not None:
            if hasattr(self, 'time'):
                time = self.time - self.time_step_output
            else:
                time = None
            if isinstance(background, xr.DataArray):
                map_x = background.coords['lon_bin']
                map_y = background.coords['lat_bin']
                scalar = background
                map_y, map_x = np.meshgrid(map_y, map_x)
            elif background == 'residence':
                scalar, lon_res, lat_res = self.get_residence_time(
                    pixelsize_m=density_pixelsize_m)
                scalar[scalar == 0] = np.nan
                lon_res, lat_res = np.meshgrid(lon_res[0:-1], lat_res[0:-1])
                lon_res = lon_res.T
                lat_res = lat_res.T
                map_x, map_y = (lon_res, lat_res)
            else:
                map_x, map_y, scalar, u_component, v_component = \
                    self.get_map_background(ax, background, crs, time=time)
                #self.time_step_output)

            if show_scalar is True:
                if contourlines is False:
                    scalar = np.ma.masked_invalid(scalar)
                    mappable = ax.pcolormesh(map_x,
                                             map_y,
                                             scalar,
                                             alpha=bgalpha,
                                             zorder=1,
                                             vmin=vmin,
                                             vmax=vmax,
                                             cmap=cmap,
                                             transform=self.crs_lonlat)
                else:
                    if contourlines is True:
                        CS = ax.contour(map_x,
                                        map_y,
                                        scalar,
                                        colors='gray',
                                        transform=self.crs_lonlat)
                    else:
                        # contourlines is an array of values
                        CS = ax.contour(map_x,
                                        map_y,
                                        scalar,
                                        contourlines,
                                        colors='gray',
                                        transform=self.crs_lonlat)
                    plt.clabel(CS, fmt='%g')

        if mappable is not None and colorbar is True:
            cb = fig.colorbar(mappable,
                              orientation='horizontal',
                              pad=cpad,
                              aspect=caspect,
                              shrink=cshrink,
                              drawedges=False)
            # TODO: need better control of colorbar content
            if clabel is not None:
                cb.set_label(clabel)
            elif isinstance(linecolor, str) and linecolor != 'gray':
                cb.set_label(str(linecolor))
            if background is not None and clabel is None:
                if isinstance(background, xr.DataArray):
                    cb.set_label(background.name)
                else:
                    cb.set_label(str(background))

        if type(background) is list:
            ax.quiver(map_x[::skip, ::skip],
                      map_y[::skip, ::skip],
                      u_component[::skip, ::skip],
                      v_component[::skip, ::skip],
                      scale=scale,
                      transform=self.crs_lonlat,
                      zorder=1)

        if lcs is not None:
            map_x_lcs, map_y_lcs = (lcs['lon'], lcs['lat'])
            ax.pcolormesh(map_x_lcs,
                          map_y_lcs,
                          lcs['ALCS'][0, :, :],
                          alpha=1,
                          vmin=vmin,
                          vmax=vmax,
                          zorder=0,
                          cmap=cmap,
                          transform=self.crs_lonlat)

        if title is not None:
            if title == 'auto':
                if hasattr(self, 'time'):
                    plt.title('%s\n%s to %s UTC (%i steps)' %
                              (self._figure_title(),
                               self.start_time.strftime('%Y-%m-%d %H:%M'),
                               self.time.strftime('%Y-%m-%d %H:%M'),
                               len(self.result.time)))
                else:
                    plt.title(
                        '%s\n%i elements seeded at %s UTC' %
                        (self._figure_title(), self.num_elements_scheduled(),
                         self.elements_scheduled_time[0].strftime(
                             '%Y-%m-%d %H:%M')))
            else:
                plt.title(title)

        if drifter is not None:
            self._plot_drifter(ax, self.crs_lonlat, drifter)

        try:
            handles, labels = ax.get_legend_handles_labels()
            if legend is not None and len(handles) > 0:
                plt.legend(loc=legend_loc, markerscale=2)
        except Exception as e:
            logger.warning('Cannot plot legend, due to bug in matplotlib:')
            logger.warning(traceback.format_exc())

        #plt.gca().tick_params(labelsize=14)
        #fig.canvas.draw()
        if filename is not None:
            plt.savefig(filename)
            logger.info('Time to make plot: ' +
                        str(datetime.now() - start_time))
        else:
            if show is True:
                plt.show()

        return ax, fig

    def _substance_name(self):
        return None

    def _figure_title(self):
        if self._substance_name() is None:
            return 'OpenDrift - ' + type(self).__name__
        else:
            return 'OpenDrift - ' + type(
                self).__name__ + ' (%s)' % self._substance_name()

    def _plot_drifter(self, ax, gcrs, drifter):
        '''Plot provided trajectory along with simulated'''
        time = np.array(drifter['time'], dtype=self.result.time.dtype)
        i = np.where((time >= self.result.time[0].values) & (time <= self.result.time[-1].values))[0]
        x, y = (np.atleast_1d(drifter['lon'])[i],
                np.atleast_1d(drifter['lat'])[i])
        dlabel = drifter['label'] if 'label' in drifter else 'Drifter'
        dcolor = drifter['color'] if 'color' in drifter else 'r'
        dlinewidth = drifter['linewidth'] if 'linewidth' in drifter else 2
        dzorder = drifter['zorder'] if 'zorder' in drifter else 10

        ax.plot(x,
                y,
                linewidth=dlinewidth,
                color=dcolor,
                transform=self.crs_lonlat,
                label=dlabel,
                zorder=dzorder)
        ax.plot(x[0], y[0], 'ok', transform=self.crs_lonlat)
        ax.plot(x[-1], y[-1], 'xk', transform=self.crs_lonlat)

    def get_map_background(self, ax, background, crs, time=None):
        # Get background field for plotting on map or animation
        # TODO: this method should be made more robust
        if time is not None:
            if isinstance(time, xr.DataArray):
                time = pd.to_datetime(time.data)
            elif isinstance(time, datetime):
                time = pd.to_datetime(time)

        if type(background) is list:
            variable = background[0]  # A vector is requested
        else:
            variable = background  # A scalar is requested
        for readerName in self.env.readers:
            reader = self.env.readers[readerName]
            if variable in reader.variables:
                if time is None or reader.start_time is None or (
                        time >= pd.Timestamp(reader.start_time)
                        and time <= pd.Timestamp(reader.end_time)) or (reader.always_valid
                                                         is True):
                    break
        if time is None:
            if hasattr(self, 'elements_scheduled_time'):
                # Using time of first seeded element
                time = self.elements_scheduled_time[0]

        # Get reader coordinates covering given map area
        axisproj = pyproj.Proj(ax.projection.proj4_params)
        xmin, xmax, ymin, ymax = ax.get_extent(self.crs_lonlat)
        cornerlons = np.array([xmin, xmin, xmax, xmax])
        cornerlats = np.array([ymin, ymax, ymin, ymax])
        reader_x, reader_y = reader.lonlat2xy(cornerlons, cornerlats)
        if sum(~np.isfinite(reader_x + reader_y)) > 0:
            # Axis corner points are not within reader domain
            reader_x = np.array([reader.xmin, reader.xmax])
            reader_y = np.array([reader.ymin, reader.ymax])
        else:
            reader_x = np.linspace(reader_x.min(), reader_x.max(), 10)
            reader_y = np.linspace(reader_y.min(), reader_y.max(), 10)

        data = reader.get_variables(background, time, reader_x, reader_y, None)
        reader_x, reader_y = np.meshgrid(data['x'], data['y'])
        if type(background) is list:  # Ensemble reader, using first member
            u_component = data[background[0]]
            v_component = data[background[1]]
            if isinstance(u_component, list):
                u_component = u_component[0]
                v_component = v_component[0]
            with np.errstate(invalid='ignore'):
                scalar = np.sqrt(u_component**2 + v_component**2)
            # NB: rotation not completed!
            u_component, v_component = reader.rotate_vectors(
                reader_x, reader_y, u_component, v_component, reader.proj,
                ccrs.PlateCarree(globe=ccrs.Globe(datum='WGS84',
                                                  ellipse='WGS84')).proj4_init)
        else:
            scalar = data[background]
            if isinstance(scalar, list):  # Ensemble reader, using first member
                scalar = scalar[0]
            u_component = v_component = None

        if reader.projected is False:
            reader_y[reader_y < 0] = 0
            reader_x[reader_x < 0] = 0

        rlons, rlats = reader.xy2lonlat(reader_x, reader_y)
        if rlons.max() > 360:
            rlons = rlons - 360
        map_x, map_y = (rlons, rlats)

        scalar = np.ma.masked_invalid(scalar)

        return map_x, map_y, scalar, u_component, v_component

    def get_lonlat_bins(self, pixelsize_m):
        latmin = self.result.lat.minval
        latmax = self.result.lat.maxval
        lonmin = self.result.lon.minval
        lonmax = self.result.lon.maxval
        deltalat = pixelsize_m / 111000.0  # m to degrees
        deltalon = deltalat / np.cos(np.radians((latmin + latmax) / 2))
        latbin = np.arange(latmin - deltalat, latmax + deltalat, deltalat)
        lonbin = np.arange(lonmin - deltalon, lonmax + deltalon, deltalon)
        return lonbin, latbin

    def get_histogram(self, pixelsize_m, **kwargs):
        from xhistogram.xarray import histogram
        lonbin, latbin = self.get_lonlat_bins(pixelsize_m)
        max_om = int(self.result.origin_marker.max().compute().values)
        origin_marker = range(max_om + 1)
        if 'weights' in kwargs and kwargs['weights'] is not None and kwargs[
                'weights'].ndim < 2:
            kwargs['weights'] = xr.DataArray(
                kwargs['weights'],
                dims=['trajectory'],
                coords={'trajectory': self.result.coords['trajectory']})
        # Xarray Dataset to store histogram per origin_marker
        h_om = xr.DataArray(np.zeros(
            (len(self.result.coords['time']), len(lonbin) - 1, len(latbin) - 1,
             max_om + 1)),
                            name='density_origin_marker',
                            dims=('time', 'lon_bin', 'lat_bin',
                                  'origin_marker'))
        h_om.coords['time'] = self.result.coords['time']
        h_om.coords['origin_marker'] = origin_marker

        for om in origin_marker:
            logger.info('\tcalculating for origin_marker %s...' % om)
            h = histogram(self.result.lon.where(self.result.origin_marker == om),
                          self.result.lat.where(self.result.origin_marker == om),
                          bins=[lonbin, latbin],
                          dim=['trajectory'],
                          **kwargs)
            if om == 0:
                h_om.coords['lon_bin'] = h.coords['lon_bin']
                h_om.coords['lat_bin'] = h.coords['lat_bin']
            h_om[:, :, :, om] = h  #.copy()
        return h_om

    def get_density_array(self, pixelsize_m, weight=None):
        lon = self.result.lon.values.T.copy()  # Copy to avoid mutating original lon/lat
        lat = self.result.lat.values.T.copy()

        deltalat = pixelsize_m / 111000.0  # m to degrees
        deltalon = deltalat / np.cos(
            np.radians((np.nanmin(lat) + np.nanmax(lat)) / 2))
        lat_array = np.arange(
            np.nanmin(lat) - deltalat,
            np.nanmax(lat) + deltalat, deltalat)
        lon_array = np.arange(
            np.nanmin(lon) - deltalat,
            np.nanmax(lon) + deltalon, deltalon)
        bins = (lon_array, lat_array)
        z = self.result.z.values.T
        if weight is not None:
            weight_array = self.result[weight].values.T

        status = self.result.status.values.T
        lon_submerged = lon.copy()
        lat_submerged = lat.copy()
        lon_stranded = lon.copy()
        lat_stranded = lat.copy()
        lon_submerged[z >= 0] = 1000
        lat_submerged[z >= 0] = 1000
        lon[z < 0] = 1000
        lat[z < 0] = 1000
        H = np.zeros((len(self.result.time), len(lon_array) - 1,
                      len(lat_array) - 1))  #.astype(int)
        H_submerged = H.copy()
        H_stranded = H.copy()
        try:
            strandnum = self.status_categories.index('stranded')
            lon_stranded[status != strandnum] = 1000
            lat_stranded[status != strandnum] = 1000
            contains_stranded = True
        except ValueError:
            contains_stranded = False

        for i in range(len(self.result.time)):
            if weight is not None:
                weights = weight_array[i, :]
            else:
                weights = None
            H[i,:,:], dummy, dummy = \
                np.histogram2d(lon[i,:], lat[i,:],
                               weights=weights, bins=bins)
            H_submerged[i,:,:], dummy, dummy = \
                np.histogram2d(lon_submerged[i,:], lat_submerged[i,:],
                               weights=weights, bins=bins)
            if contains_stranded is True:
                H_stranded[i,:,:], dummy, dummy = \
                np.histogram2d(lon_stranded[i,:], lat_stranded[i,:],
                               weights=weights, bins=bins)

        return H, H_submerged, H_stranded, lon_array, lat_array

    def get_density_array_proj(self,
                               pixelsize_m,
                               density_proj=None,
                               llcrnrlon=None,
                               llcrnrlat=None,
                               urcrnrlon=None,
                               urcrnrlat=None,
                               weight=None):
        #
        # TODO: should be merged with get_density_array
        # KFD Jan 2021
        #
        lon = self.result.lon.values.T.copy()
        lat = self.result.lat.values.T.copy()
        #deltalat = pixelsize_m/111000.0  # m to degrees
        #deltalon = deltalat/np.cos(np.radians((np.nanmin(lat) +
        #                                       np.nanmax(lat))/2))
        #lat_array = np.arange(np.nanmin(lat)-deltalat,
        #                      np.nanmax(lat)+deltalat, deltalat)
        #lon_array = np.arange(np.nanmin(lon)-deltalat,
        #                      np.nanmax(lon)+deltalon, deltalon)
        #bins=(lon_array, lat_array)
        if density_proj is None:  # add default projection with equal-area property
            density_proj = pyproj.Proj('+proj=moll +ellps=WGS84 +lon_0=0.0')
            density_proj = pyproj.Proj('+proj=longlat +a=6371229 +no_defs')

        # create a grid in the specified projection
        x, y = density_proj(lon, lat)
        if llcrnrlon is not None:
            llcrnrx, llcrnry = density_proj(llcrnrlon, llcrnrlat)
            urcrnrx, urcrnry = density_proj(urcrnrlon, urcrnrlat)
        else:
            llcrnrx, llcrnry = x.min() - pixelsize_m, y.min() - pixelsize_m
            urcrnrx, urcrnry = x.max() + pixelsize_m, y.max() + pixelsize_m

        x_array = np.arange(llcrnrx, urcrnrx, pixelsize_m)
        y_array = np.arange(llcrnry, urcrnry, pixelsize_m)
        bins = (x_array, y_array)
        outsidex, outsidey = max(x_array) * 1.5, max(y_array) * 1.5

        z = self.result.z.values.T
        if weight is not None:
            weight_array = self.result[weight].values.T

        status = self.result.status.values.T
        #lon_submerged = lon.copy()
        #lat_submerged = lat.copy()
        #lon_stranded = lon.copy()
        #lat_stranded = lat.copy()
        #lon_submerged[z>=0] = 1000
        #lat_submerged[z>=0] = 1000
        #lon[z<0] = 1000
        #lat[z<0] = 1000
        #H = np.zeros((len(times), len(lon_array) - 1,
        #              len(lat_array) - 1))#.astype(int)
        x_submerged = x.copy()
        y_submerged = y.copy()
        x_stranded = x.copy()
        y_stranded = y.copy()
        x_submerged[z >= 0] = outsidex
        y_submerged[z >= 0] = outsidey
        x[z < 0] = outsidex
        y[z < 0] = outsidey
        H = np.zeros(
            (len(self.result.time), len(x_array) - 1, len(y_array) - 1))  #.astype(int)
        H_submerged = H.copy()
        H_stranded = H.copy()
        try:
            strandnum = self.status_categories.index('stranded')
            #lon_stranded[status!=strandnum] = 1000
            #lat_stranded[status!=strandnum] = 1000
            x_stranded[status != strandnum] = outsidex
            y_stranded[status != strandnum] = outsidey
            contains_stranded = True
        except ValueError:
            contains_stranded = False

        for i in range(len(self.result.time)):
            if weight is not None:
                weights = weight_array[i, :]
            else:
                weights = None
            H[i,:,:], dummy, dummy = \
                np.histogram2d(x[i,:], y[i,:],
                               weights=weights, bins=bins)
            H_submerged[i,:,:], dummy, dummy = \
                np.histogram2d(x_submerged[i,:], y_submerged[i,:],
                               weights=weights, bins=bins)
            if contains_stranded is True:
                H_stranded[i,:,:], dummy, dummy = \
                np.histogram2d(x_stranded[i,:], y_stranded[i,:],
                               weights=weights, bins=bins)

        if density_proj is not None:
            Y, X = np.meshgrid(y_array, x_array)
            lon_array, lat_array = density_proj(X, Y, inverse=True)

        return H, H_submerged, H_stranded, lon_array, lat_array

    def get_residence_time(self, pixelsize_m):
        H,H_sub, H_str,lon_array,lat_array = \
            self.get_density_array(pixelsize_m)
        residence = np.sum(H, axis=0)
        return residence, lon_array, lat_array

    def write_netcdf_density_map(self, filename, pixelsize_m='auto'):
        '''Write netCDF file with map of particles densities'''

        if pixelsize_m == 'auto':
            lon = self.result.lon
            lat = self.result.lat
            latspan = lat.max() - lat.min()
            pixelsize_m = 30
            if latspan > .05:
                pixelsize_m = 50
            if latspan > .1:
                pixelsize_m = 300
            if latspan > .3:
                pixelsize_m = 500
            if latspan > .7:
                pixelsize_m = 1000
            if latspan > 2:
                pixelsize_m = 2000
            if latspan > 5:
                pixelsize_m = 4000

        H, H_submerged, H_stranded, lon_array, lat_array = \
            self.get_density_array(pixelsize_m)
        lon_array = (lon_array[0:-1] + lon_array[1::]) / 2
        lat_array = (lat_array[0:-1] + lat_array[1::]) / 2

        from netCDF4 import Dataset, date2num
        nc = Dataset(filename, 'w')
        nc.createDimension('lon', len(lon_array))
        nc.createDimension('lat', len(lat_array))
        nc.createDimension('time', H.shape[0])
        times = self.result.time
        if times[1] < times[0]:  # Revert for backward runs so that time is increasing
            times = times[::-1]
            H = np.flip(H, axis=0)
            H_submerged = np.flip(H_submerged, axis=0)
            H_stranded = np.flip(H_stranded, axis=0)
        timestr = 'seconds since 1970-01-01 00:00:00'
        nc.createVariable('time', 'f8', ('time', ))
        nc.variables['time'][:] = date2num(pd.to_datetime(times).to_pydatetime(), timestr)
        nc.variables['time'].units = timestr
        nc.variables['time'].standard_name = 'time'
        # Projection
        nc.createVariable('projection_lonlat', 'i8')
        nc.variables['projection_lonlat'].grid_mapping_name = \
            'latitude_longitude'
        nc.variables['projection_lonlat'].earth_radius = 6371229.
        nc.variables['projection_lonlat'].proj4 = \
            '+proj=longlat +a=6371229 +no_defs'
        # Coordinates
        nc.createVariable('lon', 'f8', ('lon', ))
        nc.createVariable('lat', 'f8', ('lat', ))
        nc.variables['lon'][:] = lon_array
        nc.variables['lon'].long_name = 'longitude'
        nc.variables['lon'].short_name = 'longitude'
        nc.variables['lon'].units = 'degrees_east'
        nc.variables['lat'][:] = lat_array
        nc.variables['lat'].long_name = 'latitude'
        nc.variables['lat'].short_name = 'latitude'
        nc.variables['lat'].units = 'degrees_north'
        # Density
        nc.createVariable('density_surface', 'u1', ('time', 'lat', 'lon'))
        H = np.swapaxes(H, 1, 2).astype('uint8')
        H = np.ma.masked_where(H == 0, H)
        nc.variables['density_surface'][:] = H
        nc.variables['density_surface'].long_name = 'Detection probability'
        nc.variables['density_surface'].grid_mapping = 'projection_lonlat'
        nc.variables['density_surface'].units = '1'
        # Density submerged
        nc.createVariable('density_submerged', 'u1', ('time', 'lat', 'lon'))
        H_sub = np.swapaxes(H_submerged, 1, 2).astype('uint8')
        H_sub = np.ma.masked_where(H_sub == 0, H_sub)
        nc.variables['density_submerged'][:] = H_sub
        nc.variables[
            'density_submerged'].long_name = 'Detection probability submerged'
        nc.variables['density_submerged'].grid_mapping = 'projection_lonlat'
        nc.variables['density_submerged'].units = '1'
        # Density stranded
        nc.createVariable('density_stranded', 'u1', ('time', 'lat', 'lon'))
        H_stranded = np.swapaxes(H_stranded, 1, 2).astype('uint8')
        H_stranded = np.ma.masked_where(H_stranded == 0, H_stranded)
        nc.variables['density_stranded'][:] = H_stranded
        nc.variables[
            'density_stranded'].long_name = 'Detection probability stranded'
        nc.variables['density_stranded'].grid_mapping = 'projection_lonlat'
        nc.variables['density_stranded'].units = '1'

        nc.close()

    def write_netcdf_density_map_proj(self,
                                      filename,
                                      pixelsize_m='auto',
                                      density_proj=None,
                                      llcrnrlon=None,
                                      llcrnrlat=None,
                                      urcrnrlon=None,
                                      urcrnrlat=None):
        '''Write netCDF file with map of particles densities for a given projection or area'''
        #
        # TODO: should be merged with write_netcdf_density_map_proj
        # KFD Jan 2021
        #

        if pixelsize_m == 'auto':
            lon = self.result.lon
            lat = self.result.lat
            latspan = lat.max() - lat.min()
            pixelsize_m = 30
            if latspan > .05:
                pixelsize_m = 50
            if latspan > .1:
                pixelsize_m = 300
            if latspan > .3:
                pixelsize_m = 500
            if latspan > .7:
                pixelsize_m = 1000
            if latspan > 2:
                pixelsize_m = 2000
            if latspan > 5:
                pixelsize_m = 4000

        if density_proj is None:  # add default projection with equal-area property
            density_proj = pyproj.Proj('+proj=moll +ellps=WGS84 +lon_0=0.0')


        H, H_submerged, H_stranded, lon_array, lat_array = \
                self.get_density_array_proj(pixelsize_m=pixelsize_m,
                                       density_proj=density_proj,
                                       llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                                       urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
        # calculate center coordinates
        logger.info(lon_array.shape, lat_array.shape)
        lon_array = (lon_array[:-1, :-1] + lon_array[1:, 1:]) / 2.
        lat_array = (lat_array[:-1, :-1] + lat_array[1:, 1:]) / 2.

        from netCDF4 import Dataset, date2num
        nc = Dataset(filename, 'w')
        nc.createDimension('x', lon_array.shape[0])
        nc.createDimension('y', lon_array.shape[1])
        nc.createDimension('time', H.shape[0])
        times = pd.to_datetime(self.result.time).to_pydatetime()
        timestr = 'seconds since 1970-01-01 00:00:00'
        nc.createVariable('time', 'f8', ('time', ))
        nc.variables['time'][:] = date2num(times, timestr)
        nc.variables['time'].units = timestr
        nc.variables['time'].standard_name = 'time'
        # Projection

        nc.createVariable('projection', 'i8')
        nc.variables['projection'].proj4 = density_proj.definition_string()
        # Coordinates
        nc.createVariable('lon', 'f8', ('y', 'x'))
        nc.createVariable('lat', 'f8', ('y', 'x'))
        nc.variables['lon'][:] = lon_array.T
        nc.variables['lon'].long_name = 'longitude'
        nc.variables['lon'].short_name = 'longitude'
        nc.variables['lon'].units = 'degrees_east'
        nc.variables['lat'][:] = lat_array.T
        nc.variables['lat'].long_name = 'latitude'
        nc.variables['lat'].short_name = 'latitude'
        nc.variables['lat'].units = 'degrees_north'

        # Density
        nc.createVariable('density_surface', 'u1', ('time', 'y', 'x'))
        H = np.swapaxes(H, 1, 2).astype('uint8')
        H = np.ma.masked_where(H == 0, H)
        nc.variables['density_surface'][:] = H
        nc.variables['density_surface'].long_name = 'Detection probability'
        nc.variables['density_surface'].grid_mapping = 'projection'
        nc.variables['density_surface'].units = '1'
        # Density submerged
        nc.createVariable('density_submerged', 'u1', ('time', 'y', 'x'))
        H_sub = np.swapaxes(H_submerged, 1, 2).astype('uint8')
        H_sub = np.ma.masked_where(H_sub == 0, H_sub)
        nc.variables['density_submerged'][:] = H_sub
        nc.variables[
            'density_submerged'].long_name = 'Detection probability submerged'
        nc.variables['density_submerged'].grid_mapping = 'projection'
        nc.variables['density_submerged'].units = '1'
        # Density stranded
        nc.createVariable('density_stranded', 'u1', ('time', 'y', 'x'))
        H_stranded = np.swapaxes(H_stranded, 1, 2).astype('uint8')
        H_stranded = np.ma.masked_where(H_stranded == 0, H_stranded)
        nc.variables['density_stranded'][:] = H_stranded
        nc.variables[
            'density_stranded'].long_name = 'Detection probability stranded'
        nc.variables['density_stranded'].grid_mapping = 'projection'
        nc.variables['density_stranded'].units = '1'

        nc.close()

    def write_geotiff(self, filename, pixelsize_km=.2):
        '''Write one GeoTiff image per timestep.

        filename should contain date identifiers, e.g. 'img_%Y%m%d_%H%M.tif'
        https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
        '''

        try:
            from osgeo import gdal, osr
        except:
            raise ValueError('GDAL is needed to write geotiff images.')
        import matplotlib.pyplot as plt
        driver = gdal.GetDriverByName('GTiff')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        colortable = gdal.ColorTable()
        colortable.SetColorEntry(0, (255, 255, 255, 0))
        colortable.SetColorEntry(1, (0, 0, 0, 255))
        colortable.SetColorEntry(2, (255, 0, 0, 255))
        colortable.SetColorEntry(3, (0, 255, 0, 255))
        colortable.SetColorEntry(4, (0, 0, 255, 255))

        lon = self.result.lon.values.T
        lat = self.result.lat.values.T
        status = self.result.status.values.T
        deltalat = pixelsize_km / 111.0  # km to degrees
        deltalon = deltalat / np.cos(np.radians((lat.min() + lat.max()) / 2))
        lat_array = np.arange(lat.min() - deltalat,
                              lat.max() + deltalat, deltalat)
        lon_array = np.arange(lon.min() - deltalat,
                              lon.max() + deltalon, deltalon)
        ilon = (np.round((lon - lon.min()) / deltalon)).astype(int)
        ilat = (np.round((lat - lat.min()) / deltalat)).astype(int)
        # Setting masked values to zero, for use as indices
        ilon[np.isinf(ilon)] = 0
        ilat[np.isinf(ilat)] = 0
        status[np.isinf(ilon)] = 0
        image = np.zeros(
            (len(self.result.time), len(lon_array), len(lat_array))).astype(int)
        geotransform = [
            lon_array.min(), deltalon, 0,
            lat_array.max(), 0, -deltalat
        ]
        for i, t in enumerate(pd.to_datetime(self.result.time)):
            image[i, ilon[i, :], ilat[i, :]] = status[i, :] + 1
            filename_i = t.strftime(filename)
            ds = driver.Create(
                filename_i,
                len(lon_array),
                len(lat_array),
                1,
                gdal.GDT_Byte,
            )
            ds.SetProjection(srs.ExportToWkt())
            ds.SetGeoTransform(geotransform)
            outband = ds.GetRasterBand(1)
            outband.SetNoDataValue(0)
            outband.WriteArray(np.fliplr(image[i, :, :]).transpose())
            outband.SetColorTable(colortable)
            ds = None

    # TODO obsolete method - to be removed
    def get_time_array(self):
        """Return a list of output times of last run."""

        logger.warning('Method get_time_array is obsolete and will soon be removed. Use o.result.time instead.')
        ## Making sure start_time is datetime, and not cftime object
        #self.start_time = datetime(self.start_time.year, self.start_time.month,
        #                           self.start_time.day, self.start_time.hour,
        #                           self.start_time.minute,
        #                           self.start_time.second)
        #td = self.time_step_output
        #time_array = [
        #    self.start_time + td * i for i in range(self.steps_output)
        #]
        #time_array_relative = [td * i for i in range(self.steps_output)]
        time_array = pd.to_datetime(self.result.time).to_pydatetime()  # Should switch to Pandas times
        time_array_relative = time_array - time_array[0]
        return time_array, time_array_relative

    def simulation_direction(self):
        """Return 1 for a forward simulation, and -1 for a backward simulation"""
        if self.time_step.days < 0:
            return -1
        else:
            return 1

    @require_mode(mode=Mode.Result)
    def plot_environment(self, filename=None, ax=None, show=True):
        """Plot mean wind and current velocities of element of last run."""
        x_wind = self.result.x_wind
        y_wind = self.result.y_wind
        wind = np.sqrt(x_wind**2 + y_wind**2)
        x_sea_water_velocity = self.result.x_sea_water_velocity
        y_sea_water_velocity = self.result.y_sea_water_velocity
        current = np.sqrt(x_sea_water_velocity**2 + y_sea_water_velocity**2)
        wind = wind.mean(dim='trajectory', skipna=True)
        current = current.mean(dim='trajectory', skipna=True)
        time = self.result.time
        time_relative = (time - time[0]) / np.timedelta64(1, 'h')

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.plot(time_relative, wind, 'b', label='Wind speed')
        ax.set_ylabel('Wind speed  [m/s]', color='b')
        ax.set_xlim([0, time_relative[-1]])
        ax.set_ylim([0, wind.max() * 1.1])

        ax2 = ax.twinx()
        ax2.plot(time_relative, current, 'r', label='Current speed')
        ax2.set_ylabel('Current speed  [m/s]', color='r')
        ax2.set_xlim([0, time_relative[-1]])
        ax2.set_ylim([0, current.max() * 1.1])
        for tl in ax.get_yticklabels():
            tl.set_color('b')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        ax.set_xlabel('Time  [hours]')
        ax.legend(loc='upper left')
        ax2.legend(loc='lower right')

        if filename is None:
            if show is True:
                plt.show()
        else:
            plt.savefig(filename)

    @require_mode(mode=Mode.Result)
    def plot_property(self, variable, filename=None, mean=False):
        """Basic function to plot time series of any output variables."""
        import matplotlib.pyplot as plt
        from matplotlib import dates

        hfmt = dates.DateFormatter('%d %b %Y %H:%M')
        fig = plt.figure()
        ax = fig.gca()
        ax.xaxis.set_major_formatter(hfmt)
        plt.xticks(rotation='vertical')

        data = self.result[variable].T
        if mean is True:  # Taking average over elements
            data = data.mean(dim='trajectory', keep_attrs=True)

        plt.plot(self.result.time, data)
        plt.title(variable)
        plt.xlabel('Time  [UTC]')
        if hasattr(data, 'units'):
            plt.ylabel('%s  [%s]' % (variable, data.units))
        else:
            plt.ylabel(variable)
        plt.subplots_adjust(bottom=.3)
        plt.grid()
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)

    @require_mode(mode=Mode.Result)
    def get_property(self, propname):
        """Get property from result, sorted by status."""
        logger.warning('Method get_property is obsolete and will soon be removed. Use o.result.<property> instead.')
        index_of_first, index_of_last = self.index_of_first_and_last()
        prop = self.result[propname].copy()
        status = self.result.status.copy()
        j = np.arange(status.shape[1])

        return prop.T, status.T

    @require_mode(mode=Mode.Result)
    def get_trajectory_lengths(self):
        """Calculate lengths and speeds along trajectories."""
        lons = self.result.lon.T
        lats = self.result.lat.T
        geod = pyproj.Geod(ellps='WGS84')
        a1, a2, distances = geod.inv(lons[0:-1, :], lats[0:-1, :],
                                     lons[1::, :], lats[1::, :])
        distances[np.isnan(distances)] = 0
        speeds = distances / self.time_step_output.total_seconds()
        distances[speeds >
                  100] = 0  # TODO: need better way to mask invalid distances
        speeds[speeds > 100] = 0  #       due to masked lons/lats arrays
        total_length = np.cumsum(distances, 0)[-1, :]

        return total_length, distances, speeds

    @require_mode(mode=Mode.Run)
    def update_positions(self, x_vel, y_vel):
        """Move particles according to given velocity components.

        This method shall account for projection metrics (a distance
        on a map projection does not necessarily correspond to the same
        distance over true ground (not yet implemented).

        Arguments:
            x_vel and v_vel: floats, velocities in m/s of particle along
                             x- and y-axes of the inherit SRS (proj4).
        """

        geod = pyproj.Geod(ellps='WGS84')

        azimuth = np.degrees(np.arctan2(x_vel, y_vel))  # Direction of motion
        velocity = np.sqrt(x_vel**2 + y_vel**2)  # Velocity in m/s
        velocity = velocity * self.elements.moving  # Do not move frosen elements

        # Calculate new positions
        self.elements.lon, self.elements.lat, back_az = geod.fwd(
            self.elements.lon, self.elements.lat, azimuth,
            velocity * self.time_step.total_seconds())

        # Check that new positions are valid
        if (self.elements.lon.min()
                < -180) or (self.elements.lon.min()
                            > 360) or (self.elements.lat.min()
                                       < -90) or (self.elements.lat.max()
                                                  > 90):
            logger.info('Invalid new coordinates:')
            logger.info(self.elements)
            sys.exit('Quitting')

    def plot_memory_usage(self, filename=None):
        plt.plot(self.memory_usage)
        plt.ylabel('Virtual memory  [GB]')
        plt.xlabel('Calculation time step')
        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)
            plt.close()

    def __repr__(self):
        """String representation providing overview of model status."""
        outStr = '===========================\n'
        if self.result is not None:
            outStr += self.performance()
            outStr += '===========================\n'
        outStr += 'Model:\t' + type(self).__name__ + \
            '     (OpenDrift version %s)\n' % opendrift.__version__
        outStr += '\t%s active %s particles  (%s deactivated, %s scheduled)\n'\
            % (self.num_elements_active(), self.ElementType.__name__,
               self.num_elements_deactivated(), self.num_elements_scheduled())
        variable_groups, reader_groups, missing = self.env.get_reader_groups()
        outStr += '-------------------\n'
        outStr += 'Environment variables:\n'
        for i, variableGroup in enumerate(variable_groups):
            outStr += '  -----\n'
            readerGroup = reader_groups[i]
            for variable in sorted(variableGroup):
                outStr += '  ' + variable + '\n'
            for i, reader in enumerate(readerGroup):
                outStr += '     ' + str(i + 1) + ') ' + reader + '\n'
        if len(self.env.missing_variables()) > 0:
            outStr += '  -----\n'
            outStr += 'Readers not added for the following variables:\n'
            for variable in sorted(self.env.missing_variables()):
                outStr += '  ' + variable + '\n'

        lazy_readers = [
            r for r in self.env.readers if self.env.readers[r].is_lazy is True
        ]
        if len(lazy_readers) > 0:
            outStr += '---\nLazy readers:\n'
            for lr in lazy_readers:
                outStr += '  ' + lr + '\n'

        outStr += '\nDiscarded readers:\n'
        for dr, reason in self.env.discarded_readers.items():
            outStr += '  %s (%s)\n' % (dr, reason)

        if hasattr(self, 'time'):
            outStr += '\nTime:\n'
            outStr += '\tStart: %s UTC\n' % (self.start_time)
            outStr += '\tPresent: %s UTC\n' % (self.time)
            if hasattr(self, 'time_step'):
                outStr += '\tCalculation steps: %i * %s - total time: %s\n' % (
                    self.steps_calculation, self.time_step,
                    self.time - self.start_time)
                outStr += '\tOutput steps: %i * %s\n' % (len(self.result.time),
                                                         self.time_step_output)
        if hasattr(self, 'messages'):
            outStr += '-------------------\n'
            outStr += self.get_messages()
        outStr += '===========================\n'
        return outStr

    def store_message(self, message):
        """Store important messages to be displayed to user at end."""
        if not hasattr(self, 'messages'):
            self.messages = []
        self.messages.append(message)

    def get_messages(self):
        """Report any messages stored during simulation."""

        if hasattr(self, 'messages'):
            return str(self.messages).strip('[]') + '\n'
        else:
            return ''

    def add_halo_readers(self):
        """Adding some Thredds and file readers in prioritised order"""

        self.add_readers_from_file(self.test_data_folder() +
                                   '../../opendrift/scripts/data_sources.txt')

    def _sphinx_gallery_filename(self, stack_offset=3):
        # This assumes that the calling script is three frames up in the stack.
        # called through a more deeply nested method stack_offset has to be changed.

        caller = inspect.stack()[stack_offset]
        caller = os.path.splitext(os.path.basename(caller.filename))[0]

        # Calling script is string input (e.g. from ..plot::)
        if caller == '<string>':
            caller = 'plot_directive'
            adir = os.path.realpath('../source/gallery/animations')
        else:
            adir = os.path.realpath('../docs/source/gallery/animations')

        if not hasattr(OpenDriftSimulation, '__anim_no__'):
            OpenDriftSimulation.__anim_no__ = {}

        if caller not in OpenDriftSimulation.__anim_no__:
            OpenDriftSimulation.__anim_no__[caller] = 0

        os.makedirs(adir, exist_ok=True)

        filename = '%s_%d.gif' % (caller,
                                  OpenDriftSimulation.__anim_no__[caller])
        OpenDriftSimulation.__anim_no__[caller] += 1

        filename = os.path.join(adir, filename)

        return filename

    def __save_animation__(self, fig, plot_timestep, filename, frames, fps,
                           blit, interval):

        if filename is None or 'sphinx_gallery' in sys.modules:
            stack_offset = 4
            import gc
            fr = sys._getframe(2)
            for o in gc.get_objects():
                if inspect.isfunction(o) and o.__code__ is fr.f_code:
                    if hasattr(getattr(self, o.__name__), '__wrapped__'):
                        stack_offset = 5
            filename = self._sphinx_gallery_filename(stack_offset=stack_offset)

        logger.info('Saving animation to ' + str(filename) + '...')

        start_time = datetime.now()

        writer = None

        if str(filename)[-4:] == '.gif':
            writer = animation.PillowWriter(fps=fps)
            # writer=animation.ImageMagickWriter(fps=fps)
        elif str(filename)[-4:] == '.mp4':
            writer = animation.FFMpegWriter(
                fps=fps,
                codec='libx264',
                bitrate=1800,
                extra_args=[
                    '-profile:v',
                    'baseline',
                    '-vf',
                    'crop=trunc(iw/2)*2:trunc(ih/2)*2',  # cropping 1 pixel if not even
                    '-pix_fmt',
                    'yuv420p',
                    '-an'
                ])
        else:
            # fallback to using funcwriter
            anim = animation.FuncAnimation(fig,
                                           plot_timestep,
                                           blit=blit,
                                           frames=frames,
                                           interval=interval)
            anim.save(filename)

        if writer is not None:
            with writer.saving(fig, filename, None):
                for i in frames if isinstance(frames, (list, range)) else range(frames):
                    plot_timestep(i)
                    writer.grab_frame()

        logger.debug(f"MPLBACKEND = {matplotlib.get_backend()}")
        logger.debug(f"DISPLAY = {os.environ.get('DISPLAY', 'None')}")
        logger.debug('Time to save animation: %s' %
                     (datetime.now() - start_time))

        plt.close()

    def calculate_ftle(self,
                       reader=None,
                       delta=None,
                       domain=None,
                       time=None,
                       time_step=None,
                       duration=None,
                       z=0,
                       RLCS=True,
                       ALCS=True):

        if reader is None:
            logger.info('No reader provided, using first available:')
            reader = list(self.env.readers.items())[0][1]
            logger.info(reader.name)
        if isinstance(reader, pyproj.Proj):
            proj = reader
        elif isinstance(reader, str):
            proj = pyproj.Proj(reader)
        else:
            proj = reader.proj

        from opendrift.models.physics_methods import ftle

        if not isinstance(duration, timedelta):
            duration = timedelta(seconds=duration)

        if domain == None:
            xs = np.arange(reader.xmin, reader.xmax, delta)
            ys = np.arange(reader.ymin, reader.ymax, delta)
        else:
            xmin, xmax, ymin, ymax = domain
            xs = np.arange(xmin, xmax, delta)
            ys = np.arange(ymin, ymax, delta)

        X, Y = np.meshgrid(xs, ys)
        lons, lats = proj(X, Y, inverse=True)

        if time is None:
            time = reader.start_time
        if not isinstance(time, list):
            time = [time]
        # dictionary to hold LCS calculation
        lcs = {'time': time, 'lon': lons, 'lat': lats}
        lcs['RLCS'] = np.zeros((len(time), len(ys), len(xs)))
        lcs['ALCS'] = np.zeros((len(time), len(ys), len(xs)))
        T = np.abs(duration.total_seconds())
        for i, t in enumerate(time):
            logger.info('Calculating LCS for ' + str(t))
            # Forwards
            if RLCS is True:
                o = self.clone()
                o.seed_elements(lons.ravel(), lats.ravel(), time=t, z=z)
                o.run(duration=duration, time_step=time_step)
                lon = o.result.lon.ffill(dim='time')
                lat = o.result.lat.ffill(dim='time')
                b_x1, b_y1 = proj(lon.T[-1].values.reshape(X.shape),
                                  lat.T[-1].values.reshape(X.shape))
                lcs['RLCS'][i, :, :] = ftle(b_x1 - X, b_y1 - Y, delta, T)
            # Backwards
            if ALCS is True:
                o = self.clone()
                o.seed_elements(lons.ravel(),
                                lats.ravel(),
                                time=t + duration,
                                z=z)
                o.run(duration=duration, time_step=-time_step)
                lon = o.result.lon.ffill(dim='time')
                lat = o.result.lat.ffill(dim='time')
                b_x1, b_y1 = proj(lon.T[-1][::-1].values.reshape(X.shape),
                                  lat.T[-1][::-1].values.reshape(X.shape))
                lcs['ALCS'][i, :, :] = ftle(b_x1 - X, b_y1 - Y, delta, T)

        lcs['RLCS'] = np.ma.masked_invalid(lcs['RLCS'])
        lcs['ALCS'] = np.ma.masked_invalid(lcs['ALCS'])
        # Flipping ALCS left-right. Not sure why this is needed
        # >> not needed anymore now that re-ordering is done above
        # lcs['ALCS'] = lcs['ALCS'][:, ::-1, ::-1]

        return lcs

    def center_of_gravity(self, onlysurface=False):
        """
        calculate center of mass and variance of all elements
        returns  (lon,lat), variance
        where (lon,lat) are the coordinates of the center of mass as
        function of time"""
        lon, lat = self.result.lon, self.result.lat
        x, y = self.proj_latlon(lon, lat)
        if onlysurface == True:
            z = self.result.z
            submerged = z < 0
            x = np.ma.array(x, mask=submerged)
            y = np.ma.array(y, mask=submerged)
        # center of gravity:
        x_m, y_m = np.ma.mean(x, axis=0), np.ma.mean(y, axis=0)
        center = self.proj_latlon(x_m, y_m, inverse=True)
        one = np.ones_like(x)
        # variance:
        variance = np.ma.mean((x - x_m * one)**2 + (y - y_m * one)**2, axis=0)

        return center, variance

    def gui_postproc(self):
        '''To be overloaded by subclasses'''
        pass
