#!/usr/bin/env python

from datetime import datetime, timedelta
from opendrift.models.leeway import Leeway
from opendrift.models.openoil3D import OpenOil3D
from opendrift.readers import reader_netCDF_CF_generic

for case in ['oil', 'leeway']:
    for timestep in [900, -900]:
        for z in [0, -200]:
            if case == 'oil':
                o = OpenOil3D(weathering_model='noaa')
                args = {'oiltype': 'HEIDRUN, STATOIL'}
            else:
                if z < 0:
                    continue  # Only surface simulation
                o = Leeway()
                args = {'objectType': 32}

            norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
            nordic = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/nordic4km/zdepths1h/aggregate_be')
            meps = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_pp_2_5km_latest.nc')
            wam = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/mywavewam4/mywavewam4_be')  
            o.add_reader([norkyst, nordic, meps, wam])

            #o.seed_elements(lon=[-0.8, 0.4], lat=[59.9, 59.95],
            #o.seed_elements(lon=[4.8, 5.1], lat=[59.9, 59.95],
            #o.seed_elements(lon=[13.11, 13.13], lat=[67.81, 67.80],
            #o.seed_elements(lon=[19.37, 19.33], lat=[70.32, 70.34],
            o.seed_elements(lon=[3.8, 3.82], lat=[59.6, 59.61], z=z,
                            time=datetime.now() - timedelta(hours=3),
                            number=1000, radius = [0, 1000], cone=True, **args)

            print o
            o.run(duration=timedelta(hours=24), time_step=timestep,
                  outfile='halo_test.nc')

            print o
            o.animation()
            if case == 'oil':
                o.plot_oil_budget()
                o.plot_property('water_fraction')
                o.animation_profile()
