#!/usr/bin/env python

from datetime import datetime, timedelta
import glob
from opendrift.models.leeway import Leeway
from opendrift.models.openoil3D import OpenOil3D
from opendrift.readers import reader_netCDF_CF_generic

readers = [  # Note that order (priority) is important!
    '/lustre/storeB/project/copernicus/sea/romsnorkyst/zdepths1h/*fc*.nc',
    'http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be',
    '/lustre/storeB/project/copernicus/sea/romsnordic/zdepths1h/roms_nordic4_ZDEPTHS_hr.fc.*.nc',
    'http://thredds.met.no/thredds/dodsC/sea/nordic4km/zdepths1h/aggregate_be',
    '/lustre/storeB/project/metproduction/products/meps/symlinks/thredds/meps_det_pp_2_5km_latest.nc',
    'http://thredds.met.no/thredds/dodsC/meps25files/meps_det_pp_2_5km_latest.nc',
    '/lustre/storeB/project/metproduction/products/arome2_5_arctic/thredds/arome_arctic_pp_2_5km_latest.nc',
    'http://thredds.met.no/thredds/dodsC/aromearcticlatest/arome_arctic_pp_2_5km_latest.nc',
    '/lustre/storeA/project/copernicus/sea/mywavewam4/*fc*.nc',
    'http://thredds.met.no/thredds/dodsC/sea/mywavewam4/mywavewam4_be',
    'http://thredds.met.no/example_of_nonexisting_url.nc',
    'example_of_nonexisting_file']

for case in ['oil', 'leeway']:  # test two models
    for timestep in [900, -900]:  # forwards and backwards
        for z in [0, -200]:  # seeding at sea surface and at 200 m depth
            if case == 'oil':
                o = OpenOil3D(weathering_model='noaa')
                args = {'oiltype': 'IVAR AASEN',
                #args = {'oiltype': 'WISTING',
                        'z': z}
            else:
                if z < 0:
                    continue  # Only surface seeding for Leeway
                o = Leeway()
                args = {'objectType': 32}

            o.add_readers_from_list(readers, timeout=5)
            print(o)

            #lons=[-0.8, 0.4]; lats=[59.9, 59.95] # NorKyst border
            #lons=[4.8, 5.1]; lats=[59.9, 59.95] # Western Norway coast
            #lons=[13.11, 13.13]; lats=[67.81, 67.80] # Lofoten
            #lons=[19.37, 19.33]; lats=[70.32, 70.34] # Troms
            lons=[3.8, 3.82]; lats=[59.6, 59.61] # North Sea

            o.seed_elements(lon=lons, lat=lats,
                            time=[datetime.now() - timedelta(hours=3),
                                  datetime.now()],
                            number=1000, radius = [0, 1000],
                            cone=True, **args)

            print(o)
            o.run(duration=timedelta(hours=24), time_step=timestep,
                  time_step_output=1800, outfile='halo_test.nc')

            print(o)
            o.animation()
            if case == 'oil':
                o.plot()
                o.plot_oil_budget()
                o.plot_property('water_fraction')
                o.animation_profile()
