#!/usr/bin/env python
"""
Using wind measurements
=======================
"""

from datetime import datetime, timedelta, date
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.openoil import OceanDrift
from opendrift.readers import reader_timeseries
from opendrift.readers.operators import ops
import opendrift.models.physics_methods as phy
import xarray as xr
import trajan as _
import numpy as np
import pyproj
import matplotlib.pyplot as plt
import matplotlib as mpl

# Preparation of figsize, dpi, fontsize and font style
# Could also be done using mpl style sheet
# https://matplotlib.org/stable/users/explain/customizing.html
mpl.rcParams['figure.figsize'] = (5, 5)
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "DejaVu Sans"
mpl.rcParams['font.size'] = 22
mpl.rcParams['mathtext.fontset'] = "dejavusans"



#%%
# This example shows how to use wind or currents measurments 
# combined with an ocean model to improove the drift
# calculation arround the measurment station.


#%% 
# Fisrt, we import two times the same model, because one
# will be used for the usual simulation, and one will
# be used for the measurement combined simulation.

reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
reader_norkyst2 = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

#Verification
print(reader_norkyst)


#%% 
# We then import wind measurement data,
# that will be stored in a timeserie reader.
# Don't forget to give the lon and lat of the 
# measurement station.

time_list = [datetime(2024, 5, 27, 6, 00, 00) + i*timedelta(minutes=10) for i in range(7)] # Manual import for this example
wind_speed_array = np.array([1.9, 3.0, 3.5, 3.4, 3.9, 2.8, 2.5])
wind_speed_angle = np.array([22, 25, 24, 19, 17, 30, 16])



# Alternative if data are in a csv file

# data_wind_bridge = np.genfromtxt("Simulation_data/data_wind_21_mai_2024.txt", dtype = '<U10', skip_header=1)
# temporal_time = np.array([data_wind_bridge[i, 0] + "T" + data_wind_bridge[i, 1] for i in range(len(data_wind_bridge))])
# time_list = [datetime.fromisoformat(temporal_time[i]) for i in range(len(temporal_time))]
# wind_speed_array = data_wind_bridge[:, 3]
# wind_speed_angle = data_wind_bridge[:, 5]
# wind_speed_array = wind_speed_array.astype(float)
# wind_speed_angle = wind_speed_angle.astype(float)


#Calculation of wind_x and wind_y
wind_speed_angle = np.deg2rad(wind_speed_angle)
wind_x = wind_speed_array * np.cos(wind_speed_angle)
wind_y = wind_speed_array * np.sin(wind_speed_angle)

#Creation of the wind measurement reader
lon, lat = 5.17, 60.37
parameter_value_map = {'time': time_list, 'wind_speed': wind_speed_array, 'x_wind': wind_x, 'y_wind':  wind_y, 'lon':lon, 'lat':lat}
wind_obs_reader = reader_timeseries.Reader(parameter_value_map)

#Preparing future plots
text = [{'s': '* Station', 'x': lon, 'y': lat, 'fontsize': 20, 'color': 'k', 'backgroundcolor': 'white', 'bbox': dict(facecolor='none', edgecolor='none', alpha=0.4), 'zorder': 1000}]

#Verification
print(wind_obs_reader)

#%% 
# Now we combine the second norkyst reader 
# with the measurement one, using a gaussian 
# combined reader.

c = reader_norkyst2.combine_gaussian(wind_obs_reader, std = 1000)
print(c)


#%%
# Configuration of the simulations

time_start = datetime(2024, 5, 27, 6, 00, 00)
time_end = time_start + timedelta(hours = 2) 
time = [time_start, time_end]
time_step = 150
time_step_output = 150

# Seeding configuration
lon = 5.187
lat = 60.385
radius = 0
number = 1
wind_drift_factor = 0.03



#%% 
# Usual simulation

o = OceanDrift(loglevel=20)
o.add_reader([reader_norkyst])

# o.set_config('drift:current_uncertainty', .1)
# o.set_config('drift:wind_uncertainty', 1)

o.seed_elements(lon=lon, lat=lat, radius=radius, number=number, time=time_start, wind_drift_factor=wind_drift_factor)
o.run(end_time=time_end, time_step=time_step, time_step_output=time_step_output)
o.plot(fast = False, background=['x_sea_water_velocity', 'y_sea_water_velocity'], legend=['Norkyst only', 'Gaussian measurement'], buffer = .023, markersize = 150, linewidth = 3, title = "", xlocs = mpl.ticker.MaxNLocator(5), ylocs = mpl.ticker.MaxNLocator(5), clabel = r"Wind speed $\mathrm{(m.s^{-1})}$", cpad = 0.08, text=text)

#%%
# Measurement modified model simulation
# and comparison
o2 = OceanDrift(loglevel=20)
o2.set_config('drift:max_speed', 10) #Necessary for the moment because of a bug to be solved
o2.add_reader([c, reader_norkyst])

#o2.set_config('drift:current_uncertainty', .1)
#o2.set_config('drift:wind_uncertainty', 1)

o2.seed_elements(lon=lon, lat=lat, radius=radius, number=number, time=time_start, wind_drift_factor=wind_drift_factor)
o2.run(end_time=time_end, time_step=time_step, time_step_output=time_step_output)

o.plot(fast = False, compare=o2, background=['x_sea_water_velocity', 'y_sea_water_velocity'], legend=['Norkyst only', 'Gaussian measurement'], buffer = .023, markersize = 150, linewidth = 3, title = "", xlocs = mpl.ticker.MaxNLocator(5), ylocs = mpl.ticker.MaxNLocator(5), clabel = r"Wind speed $\mathrm{(m.s^{-1})}$", cpad = 0.08, text=text)

#%%
# Here, we generate more particles and look
# at the skillscore impact of the measurment 
# combined simulation
    
radius = 100
number = 100
  
o = OceanDrift(loglevel=20)
o.add_reader([reader_norkyst])
o.seed_elements(lon=lon, lat=lat, radius=radius, number=number, time=time_start, wind_drift_factor=wind_drift_factor)
o.run(end_time=time_end, time_step=time_step, time_step_output=time_step_output)

o2 = OceanDrift(loglevel=20)
o2.set_config('drift:max_speed', 10) #Necessary for the moment because of a bug to be solved
o2.add_reader([c, reader_norkyst])
o2.seed_elements(lon=lon, lat=lat, radius=radius, number=number, time=time_start, wind_drift_factor=wind_drift_factor)
o2.run(end_time=time_end, time_step=time_step, time_step_output=time_step_output)

o.plot(fast = False, compare=o2, background=['x_sea_water_velocity', 'y_sea_water_velocity'], legend=['Norkyst only', 'Gaussian measurement'], buffer = .023, markersize = 70, linewidth = 1, title = "", xlocs = mpl.ticker.MaxNLocator(5), ylocs = mpl.ticker.MaxNLocator(5), clabel = r"Wind speed $\mathrm{(m.s^{-1})}$", cpad = 0.08, text=text)

#%%
# Calculating skillscores with TrajAn
skillscores = o.result.traj.skill(o2.result, method='liu-weissberg', tolerance_threshold=1)

plt.figure()
plt.hist(skillscores, bins = 100, range = (0, 1), facecolor = 'none', edgecolor = 'C0')
plt.xlabel("Skillscore")
plt.ylabel("Histogram (percentages)")
plt.show()

plt.figure()
plt.hist(skillscores, bins = 100, range = (0.9, 0.95), facecolor = 'none', edgecolor = 'C0')
plt.xlabel("Skillscore")
plt.ylabel("Histogram (percentages)")
plt.show()
