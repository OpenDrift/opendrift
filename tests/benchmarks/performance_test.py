#!/usr/bin/env python

import os
import sys
from datetime import datetime
try:
    from psutil import virtual_memory
    ram = virtual_memory().total/(1024**3)
except:
    ram = 'unknown'
import multiprocessing
import numpy as np
import matplotlib
import scipy
import platform
import netCDF4
from opendrift.models.openoil import OpenOil
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_global_landmask
from opendrift.readers.interpolation import ReaderBlock

print('\n\n')
print('==============================================')
print('=         OpenDrift performance test         =')
print('==============================================')
print('Reference machine:')
print('  Ubuntu PC with 8 GB memory and')
print('  8 processors (Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz)')
print('  NumPy version 1.11.2')
print('  SciPy version 0.9.0')
print('  Matplotlib version 2.0.0')
print('  NetCDF4 version 1.0.7')
print('  Python version 2.7.3 (default, Jun 22 2015, 19:33:41) [GCC 4.6.3]')
print('  Graphics card:')
print('    01:00.0 VGA compatible controller: NVIDIA Corporation GT218 [NVS 300] (rev a2)')

print('------------------------------------------------')
print('This machine:')
print('  %s GB memory' % ram)
print('  %s processors (%s)' % (multiprocessing.cpu_count(),
                                platform.processor()))
print('  NumPy version %s' % np.__version__)
print('  SciPy version %s' % scipy.__version__)
print('  Matplotlib version %s' % matplotlib.__version__)
print('  NetCDF4 version %s' % netCDF4.__version__)
print('  Python version %s' % sys.version.replace('\n', ''))
print('  Graphics card:')
try:
    print('    ' + os.popen('lspci | grep VGA').read().strip())
except:
    print('    unknown')
print('------------------------------------------------')
print('\n')

print('Test 1: generation of landmask instance at full resolution')
print('  54.0 seconds on reference machine.')
start_time = datetime.now()
reader_landmask = reader_global_landmask.Reader()
time_spent = datetime.now() - start_time
print('%6.1f seconds on this machine' % time_spent.total_seconds())


print('------------------------------------------------')
print('Test 2: Reading from netCDF file')
print('  0.05 seconds on reference machine.')
o = OpenOil(loglevel=50) # Quiet
reader_arctic = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '2Feb2016_Nordic_sigma_3d/Arctic20_1to5Feb_2016.nc')
x = reader_arctic.x[10:12]
y = reader_arctic.y[10:12]
z = np.array([-20, -10])
variables = ['sea_surface_elevation', 'sea_ice_area_fraction', 'y_sea_water_velocity', 'sea_floor_depth_below_sea_level', 'sea_water_salinity', 'x_sea_water_velocity', 'sea_water_temperature', 'sea_ice_thickness']
reader_arctic.buffer=1000 # read all
reader_arctic.verticalbuffer=1000 # read all
start_time = datetime.now()
for t in reader_arctic.times:
    v = reader_arctic.get_variables(variables, time=t, x=x, y=y, z=z)
time_spent = datetime.now() - start_time
print('%6.2f seconds on this machine' % time_spent.total_seconds())


print('--------------------------------------------------------')
print('Test 3: Interpolating 3D arrays onto particle positions')
print('   2.5 seconds on reference machine.')
start_time = datetime.now()
b = ReaderBlock(v, interpolation_horizontal='linearND')
num_points = 10000
x = np.random.uniform(reader_arctic.xmin, reader_arctic.xmax, num_points)
y = np.random.uniform(reader_arctic.ymin, reader_arctic.ymax, num_points)
z = np.random.uniform(-200, 0, num_points)
env, prof = b.interpolate(x, y, z, variables,
                          profiles=['sea_water_temperature'],
                          profiles_depth=[-30, 0])
time_spent = datetime.now() - start_time
print('%6.1f seconds on this machine' % time_spent.total_seconds())


print('--------------------------------------------------------')
print('Test 4: Vertical mixing with 50 elements and 7200 cycles (CPU-heavy)')
print('  10.0 seconds on reference machine.')
reader_arctic.buffer=10
reader_arctic.verticalbuffer=1
o = OpenOil(loglevel=50) # Quiet
o.add_reader(reader_arctic)
o.set_config('environment:fallback:x_wind', 10)
o.set_config('environment:fallback:y_wind', 0)
o.set_config('vertical_mixing:timestep', 1)
o.seed_elements(lon=15, lat=72, number=50, radius=10000,
                time=reader_arctic.start_time)
start_time = datetime.now()
o.run(steps=2)
time_spent = datetime.now() - start_time
print('%6.1f seconds on this machine' % time_spent.total_seconds())


print('--------------------------------------------------------')
print('Test 5: Vertical mixing with 500000 elements and 10 cycles (memory-heavy)')
print('  38.0 seconds on reference machine.')
o = OpenOil(loglevel=50) # Quiet
o.add_reader(reader_arctic)
o.set_config('environment:fallback:x_wind', 10)
o.set_config('environment:fallback:y_wind', 0)
o.set_config('vertical_mixing:timestep', 50)
o.seed_elements(lon=15, lat=72, number=500000, radius=10000,
                time=reader_arctic.start_time)
start_time = datetime.now()
o.run(steps=1, time_step=500)
time_spent = datetime.now() - start_time
print('%6.1f seconds on this machine' % time_spent.total_seconds())



print('\n\n')
