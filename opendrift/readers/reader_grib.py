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

from datetime import datetime
import logging
logger = logging.getLogger(__name__)

import numpy as np
try:
    import pygrib
except:
    raise ImportError('PyGrib library is needed for GRIB files: '
                      'https://jswhit.github.io/pygrib/docs/index.html')

from opendrift.readers.basereader import BaseReader, StructuredReader

# Hardcoded "GRIB-tables" for now.
grib_variable_mapping = {
    'enmi': {  # MET Norway (enmi) GRIB codes (marsParam)
        '33.3': 'x_wind',
        '34.3': 'y_wind',
        '49.3': 'x_sea_water_velocity',
        '50.3': 'y_sea_water_velocity',
        '100.3': 'sea_surface_wave_significant_height',
        '231.3':
            'sea_surface_wave_period_at_variance_spectral_density_maximum',
        '232.3':
            'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment',
        '247.3': 'sea_surface_wave_stokes_drift_x_velocity',
        '248.3': 'sea_surface_wave_stokes_drift_y_velocity'},
    'ecmf': {  # ECMWF GRIB codes
        '165.128': 'x_wind',
        '166.128': 'y_wind',
        '215.140': 'sea_surface_wave_stokes_drift_x_velocity',
        '216.140': 'sea_surface_wave_stokes_drift_y_velocity',
        '229.140': 'sea_surface_wave_significant_height',
        '231.140':
            'sea_surface_wave_period_at_variance_spectral_density_maximum',
        '232.140':
            'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment'
        },
    'kwbc': {  #
        '33.2': 'x_wind',
        '34.2': 'y_wind',
        '49.2': 'x_sea_water_velocity',
        '50.2': 'y_sea_water_velocity'}
     }


class Reader(BaseReader, StructuredReader):

    def __init__(self, filename=None, name=None):

        if filename is None:
            raise ValueError('Need filename as argument to constructor')

        if name is None:
            self.name = filename
        else:
            self.name = name

        try:
            # Open file, check that everything is ok
            logger.info('Opening dataset: ' + filename)
            self.grib = pygrib.open(filename)
        except:
            raise ValueError('Could not open ' + filename +
                             ' with pygrib library')

        ####################
        # Scan file
        ####################
        levels = []
        centre = []
        projs = []
        times = []
        marsParams = []
        for m in self.grib:
            levels.append(m.level)
            centre.append(m.centre)
            projs.append(m.projparams['proj'])
            s = '%s%04d' % (m.validityDate, m.validityTime)
            times.append(datetime.strptime(s, '%Y%m%d%H%M'))
            marsParams.append(m.marsParam)

        ################
        # Projection
        ################
        projs = list(set(projs))
        if len(projs) > 1:
            raise ValueError('File with data in several projections is not '
                             'supported: ' + str(projs))
        if projs[0] == 'cyl':
            self.proj4 = '+proj=latlong'
        else:
            raise ValueError('Only GRIB files with latlon-projection are '
                             'currently supported, given projection is: %s' %
                             m.projparams)
        lalo = m.latlons()
        x = m.distinctLongitudes
        y = m.distinctLatitudes
        self.xmin = x.min()
        self.xmax = x.max()
        self.ymin = y.min()
        self.ymax = y.max()
        self.delta_x = x[1] - x[0]
        self.delta_y = y[1] - y[0]

        ####################################
        # GRIB source and parameter names
        ####################################
        centre = list(set(centre))
        if len(centre) > 1:
            raise ValueError('File contains data from several centres: ' +
                             str(centre))
        else:
            centre = centre[0]
        if centre in grib_variable_mapping:
            self.grib_mapping = grib_variable_mapping[centre]
            logger.info('Parameter codes not defined in mapper: ' +
                         str(set(marsParams) - set(self.grib_mapping)))
        else:
            raise ValueError(
                'No GRIB variable mapping defined for centre ' + centre)
        self.marsParams = list(set(self.grib_mapping) &
                               set(marsParams))
        self.variables = [self.grib_mapping[v] for v in self.marsParams]

        ####################################
        # Sort by variable and time
        ####################################
        levels = np.array(levels)
        self.indices = {}
        self.levels = {}
        for i, var in enumerate(self.variables):
            m = list(self.grib_mapping.keys())[list(self.grib_mapping.values()).index(var)]
            self.indices[var] = np.where(np.array(marsParams) == m)[0]
            self.levels[var] = levels[self.indices[var]]
        if len(self.variables) > 0:
            self.times = [times[k] for k in self.indices[var]]
            if len(self.times) > 1 and self.times[0] == self.times[1]:  # Duplicate times
                logger.info('Duplicate times for variables, using '
                             'only first occurance.')
                lasttime = None
                for i, tim in enumerate(self.times):
                    if 'newind' not in locals():
                        newind = [i]
                        lasttime = tim
                    else:
                        if tim != lasttime:
                            newind.append(i)
                            lasttime = tim
                newind = np.array(newind, dtype=np.int32)
                self.times = [self.times[i] for i in newind]
                for i, var in enumerate(self.variables):  # Updating indices
                    self.indices[var] = self.indices[var][newind]
                    self.levels[var] = levels[self.indices[var]]
            self.start_time = self.times[0]
            self.end_time = self.times[-1]
            if len(self.times) > 1:
                self.time_step = self.times[-1] - self.times[-2]
            else:
                self.time_step = None
        else:
            self.start_time = None
            self.end_time = None
            self.time_step = None

        # Run constructor of parent Reader class
        super(Reader, self).__init__()

    def get_variables(self, requested_variables, time=None,
                      x=None, y=None, z=None):

        requested_variables, time, x, y, z, outside = self.check_arguments(
            requested_variables, time, x, y, z)

        nearestTime, dummy1, dummy2, indxTime, dummy3, dummy4 = \
            self.nearest_time(time)

        variables = {}
        delta = self.buffer*self.delta_x
        lonmin = np.maximum(x.min() - delta, self.xmin)
        lonmax = np.minimum(x.max() + delta, self.xmax)
        latmin = np.maximum(y.min() - delta, self.ymin)
        latmax = np.minimum(y.max() + delta, self.ymax)

        for var in requested_variables:
            ind = int(self.indices[var][indxTime]) + 1
            msg = self.grib[ind]
            variables[var], lats, lons = msg.data(lat1=latmin, lat2=latmax,
                                                  lon1=lonmin, lon2=lonmax)
        variables['x'] = lons[0, :]
        variables['y'] = lats[:, 0]
        variables['z'] = None
        variables['time'] = nearestTime

        return variables
