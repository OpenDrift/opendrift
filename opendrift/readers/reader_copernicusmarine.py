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
# Copyright 2024, Knut-Frode Dagestad, MET Norway

import logging
logger = logging.getLogger(__name__)
from opendrift.readers.reader_netCDF_CF_generic import Reader as Reader_CF_generic

class Reader(Reader_CF_generic):
    '''A wrapper around reader_netCDF_CF_generic for CMEMS datasets'''

    def __init__(self, dataset_id, username=None, password=None):

        try:
            import copernicusmarine
        except:
            raise ValueError('Copernicus Marine Client is not installed')

        if username is None:
            try:
                import os
                import netrc
                try:
                    n = netrc.netrc()
                    netrcfile = os.path.join(os.path.expanduser('~'), '.netrc')
                except Exception as e:
                    import opendrift
                    netrcfile = os.path.join(os.path.dirname(opendrift.__file__), '.netrc')
                    n = netrc.netrc(os.path.join(netrcfile))
                machines = ['copernicusmarine', 'nrt.cmems-du.eu']
                for m in machines:
                    try:
                        username, dummy, password = n.authenticators(m)
                        logger.info(f'Using CMEMS password for user {username} from {netrcfile} [{m}]')
                        break
                    except:
                        pass
            except:
                pass

        if username is None:
            raise ValueError('To use CMEMS datasets, provide username and password, or store these in .netrc file in home folder or main opendrift folder with machine name "copernicusmarine"')

        ds = copernicusmarine.open_dataset(dataset_id=dataset_id, username=username, password=password)
        ds['name'] = str(dataset_id)

        # Run constructor of parent Reader class
        super(Reader, self).__init__(ds, name=dataset_id)
