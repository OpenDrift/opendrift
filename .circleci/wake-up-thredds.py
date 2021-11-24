#!/usr/bin/env python

import urllib.request as urllib_request
import time
from netCDF4 import Dataset

TIMEOUT = 5.
RETRY_DELAY = 2.
TRIES = 3

for r in [
        'https://thredds.met.no/thredds/users.html',
        'https://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/netcdf_unstructured_samples/AkvaplanNiva_sample_lonlat_fixed.nc.html',
        'https://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/netcdf_unstructured_samples/schism_marl20080101_00z_3D.nc.html',
        'http://tds0.ifremer.fr',
        'http://tds0.ifremer.fr/thredds/dodsC/GLOBCURRENT-L4-CUREUL_15M-ALT_SUM_NRT-V03.0.html',
        'https://iws.ismar.cnr.it/thredds/dodsC/emerge/shyfem_unstructured_adriatic.nc.html'
]:
    print('Waking up: %s..' % r)
    for _ in range(TRIES):
        try:
            # Wake up thredds-server, if sleeping
            request = urllib_request.Request(r)
            assert urllib_request.urlopen(request, timeout=TIMEOUT).getcode() == 200
            print("ok!")

            if '.nc.html' in r:
                ncurl = r[:-5]
                print('testing DAP: %s..' % ncurl)
                d = Dataset(ncurl)
                assert len(d.variables) > 0
                print('ok!')

            break
        except:
            print("failed, trying again..")
            time.sleep(RETRY_DELAY)
            pass
        print("failed, giving up.")

