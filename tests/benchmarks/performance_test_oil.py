#!/usr/bin/env python

import sys
from datetime import datetime, timedelta
from opendrift.models.openoil3D import OpenOil3D
from opendrift.readers import reader_netCDF_CF_generic

if __name__ == '__main__':

    if len(sys.argv) == 1:
        mode = 'fast'
    if len(sys.argv) == 2:
        mode = sys.argv[1]
        if mode not in ['fast', 'slow']:
            sys.exit('Usage: %s [fast|slow]' % sys.argv[0])
    if len(sys.argv) > 2:
        sys.exit('Usage: %s [fast|slow]' % sys.argv[0])

    o = OpenOil3D(weathering_model='default')

    readers_thredds = [  # Note that order (priority) is important!
        'http://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/case_1_mar_2017/NorKyst-800m_ZDEPTHS_his_00_1March2017.nc',
        'http://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/case_1_mar_2017/Nordic-4km_ZDEPTHS_avg_00.nc',
        'http://thredds.met.no/thredds/dodsC/metusers/knutfd/thredds/case_1_mar_2017/meps_mbr0_pp_2_5km_20170301T00Z.nc',
        'http://thredds.met.no//thredds/dodsC/metusers/knutfd/thredds/case_1_mar_2017/MyWave_wam4_WAVE_20170301T00Z.nc',
        'http://thredds.met.no//thredds/dodsC/metusers/knutfd/thredds/case_1_mar_2017/MyWave_wam4_RAD_20170301T00Z.nc',
        ]

    readers_local = [  # Note that order (priority) is important!
        o.test_data_folder() + '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc',
        o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc']


    if mode == 'slow':
        readers = readers_thredds
        start_time = datetime(2017, 3, 1, 0)
        number = 20000
        hours = 60
        source = 'thredds'
        o.set_config('general:basemap_resolution', 'f')
    elif mode == 'fast':
        readers = readers_local
        start_time = datetime(2015, 11, 16, 0)
        number = 10000
        hours = 24
        source = 'file'
        o.set_config('general:basemap_resolution', 'h')

    print('#'*40)
    print('Running %s test:\n %s elements\n %s hours\n input data from %s' % (
                mode, number, hours, source))
    print('#'*40)

    o.add_readers_from_list(readers, timeout=5)
    o.seed_elements(lon=4.6, lat=59.7, time=start_time,
                    number=number, radius=3000)

    print(o)
    o.run(duration=timedelta(hours=hours), time_step=900,
          time_step_output=3600, outfile='performance_test_oil.nc')

    print(o)
    #o.plot()
    #o.animation()
    #o.animation_profile()
