#!/usr/bin/env python
"""
CMEMS current components
========================

The global CMEMS/copernicus ocean model contains several different surface current fields.
OpenDrift will by default use the variables with "standard_name" attribute equal to what is requested from a module, typically "x_sea_water_velocity"/"y_sea_water_velocity" for currents.
Note that "east/north" counterparts will also be detected, and eventual rotation will be performed automatically.

This example illustrates how a "standard_name_mapping" can be added to the generic netCDF reader to chose alternative variables.
The example also illustrates the alternative (experimental) mechanism of summing two readers.
"""

from datetime import datetime, timedelta
import cf_xarray
import copernicusmarine
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers.reader_netCDF_CF_generic import Reader

#%%
# Get an Xarray dataset from copernicusmarine client
ds = copernicusmarine.open_dataset(dataset_id='cmems_mod_glo_phy_anfc_merged-uv_PT1H-i', chunk_size_limit=0)
print(ds)     # Default Xarray output
print(ds.cf)  # Output from cf-xarray

#%%
# Create an OpenDrift reader from this dataset
reader_default = Reader(ds, name='CMEMS default')

#%%
# Mapping other variables to required standard_name's
reader_tides = Reader(ds, standard_name_mapping={
                        'utide': 'x_sea_water_velocity',
                        'vtide': 'y_sea_water_velocity',
                        }, name='Tides only')
reader_stokes = Reader(ds, standard_name_mapping={
                        'vsdx': 'x_sea_water_velocity',
                        'vsdy': 'y_sea_water_velocity',
                        }, name='Stokes only')
reader_total = Reader(ds, standard_name_mapping={
                        'utotal': 'x_sea_water_velocity',
                        'vtotal': 'y_sea_water_velocity',
                        }, name='Total current')

#%%
# Run and compare simulations using the different current components
cases = {'Eulerian current': reader_default,
         'Tides only': reader_tides,
         'Stokes drift only': reader_stokes,
         'Total current': reader_total,
         'SUM: Eulerian + Tides': reader_default + reader_tides,  # Experimental feature
         'SUM: Eulerian + Stokes': reader_default + reader_stokes,     # Experimental feature
         'SUM: Eulerian + Tides + Stokes': reader_default + reader_tides + reader_stokes  # Experimental feature
         }

simulations = []
for cname, reader in cases.items():
    o = OceanDrift()
    o.add_reader(reader, variables=['x_sea_water_velocity', 'y_sea_water_velocity'])
    o.seed_elements(lon=4.8, lat=60, time=datetime(2024, 10, 31, 6))
    o.run(duration=timedelta(hours=12))
    simulations.append(o)

simulations[0].plot(filename='cmems_comparison.png', buffer=.05,
                    compare=simulations[1:], legend=list(cases))
