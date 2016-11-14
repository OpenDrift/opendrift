#!/usr/bin/env python

from datetime import datetime, timedelta
from opendrift.models.leeway import Leeway
from opendrift.models.openoil3D import OpenOil3D
from opendrift.readers import reader_netCDF_CF_generic

b = 'h'
for case in ['oil', 'leeway']:
    for timestep in [900, -900]:
        if case == 'oil':
            o = OpenOil3D(weathering_model='noaa', basemap_resolution=b)
            args = {'oiltype': 'HEIDRUN, STATOIL'}
        else:
            o = Leeway(basemap_resolution=b)
            args = {'objectType': 32}

        norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
        nordic = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/nordic4km/zdepths1h/aggregate_be')
        arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
        #norkyst.interpolation = 'linearND'
        o.add_reader([norkyst, nordic, arome])
        #o.add_reader([nordic, arome])

        #o.seed_elements(lon=[-0.8, 0.4], lat=[59.9, 59.95],
        #o.seed_elements(lon=[4.8, 5.1], lat=[59.9, 59.95],
        #o.seed_elements(lon=[13.11, 13.13], lat=[67.81, 67.80],
        #o.seed_elements(lon=[19.37, 19.33], lat=[70.32, 70.34],
        o.seed_elements(lon=[3.8, 3.82], lat=[59.6, 59.61],
                        time=datetime.now() - timedelta(hours=3),
                        number=1000, radius = [0, 1000], cone=True, **args)

        print o
        o.run(duration=timedelta(hours=48), time_step=timestep,
              outfile='halo_test.nc')

        print o
        o.animation()
        if case == 'oil':
            o.plot_oil_budget()
            o.animation_profile()
