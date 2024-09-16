#!/usr/bin/env python
"""
Icebergs (openberg)
====================
************************ Data sources *********************
*. TOPAZ5: 2021-07-05  to  2024-03-21   cmems_mod_arc_phy_anfc_6km_detided_P1D-m   ###Can be useful  cmems_mod_arc_phy_anfc_6km_detided_PT1H-i
*. TOPAZ6: 2018-01-01  to  2024-03-21   dataset-topaz6-arc-15min-3km-be   ###2D tides model
*. TOPAZ4: 1991-01-01  to  2022-12-31   cmems_mod_arc_phy_my_topaz4_P1D-m   ###No more used
*. GLORYS: 1993-01-01  to  2022-07-31   cmems_mod_glo_phy_myint_0.083deg_P1D-m  or  cmems_mod_glo_phy_my_0.083deg_P1D-m   ###Can be useful
*. Globcurrent: 1993-01-01  to  2022-07-31   cmems_obs_mob_glo_phy-cur_my_0.25deg_P1D-m   ###Sat. product
*. Waves: 1977-01-01  to  2022-11-30 cmems_mod_arc_wav_my_3km_PT1H-i  
*. Waves: start: 2021-10-01 00:00:00   end: 2024-09-21 12:00:00  https://thredds.met.no/thredds/dodsC/cmems/mywavewam3km/dataset-wam-arctic-1hr3km-be.ncml
"""
#import copernicus_marine_client as copernicusmarine
from opendrift.readers.reader_netCDF_CF_generic import Reader
from opendrift.models.openberg import OpenBerg
from datetime import datetime,timedelta
import numpy as np
import os
import matplotlib.pyplot as plt
import opendrift



# # Credentials of CMEMS Account
# USERNAME='USERNAME'
# PASSWORD='PASSWORD'

# # Define datasets
# Wind = '/Volumes/ACH_HDD/DATA_ACCIBERG/Wind/6h.10m_wind_2021.nc'
# Ocean = 'cmems_mod_glo_phy_myint_0.083deg_P1D-m'
# Wave = 'cmems_mod_glo_wav_my_0.2deg_PT3H-i'
# # Add Readers
# reader_wind = Reader(Wind)
# ##print(reader_wind)

# ds_current = copernicusmarine.open_dataset(dataset_id=Ocean, username=USERNAME, password=PASSWORD)
# reader_current = Reader(ds_current)
# ##print(reader_current)

# ds_wave = copernicusmarine.open_dataset(dataset_id=Wave, username=USERNAME, password=PASSWORD)
# reader_wave = Reader(ds_wave)
# ##print(reader_wave)
'''
o = OpenBerg()
o.set_config('drift:wave_rad', False)
o.set_config('drift:stokes_drift', False)
o.set_config('drift:coriolis', False)
o.set_config('drift:vertical_profile', True)
o.set_config('processes:grounding', False)
o.set_config('processes:roll_over', False)
#o.set_config('processes:melting', False)

##### Dataset from copernicus_marine_client
#import copernicus_marine_client as copernicusmarine
#ds = copernicusmarine.open_dataset(
#    dataset_id='cmems_obs_mob_glo_phy-cur_my_0.25deg_P1D-m', 
#    username='aothmani', password='ach@REF4098')
#r = Reader(ds)
#o.add_reader(r)



#o.add_readers_from_list([
#        'https://thredds.met.no/thredds/dodsC/cmems/topaz6/dataset-topaz6-arc-15min-3km-be.ncml',
#        'https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd'])

reader_current = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/TP5Metno_m09_m10_y2021/2021*_dm-metno-MODEL-topaz5-ARC-b2021*-fv02.0.nc') ###TP5 Metno 3D
#reader_wind = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/Wind/6h.10m_wind_2021.nc')
#reader_wave = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/Waves/MyWam3km_hindcast-cmems_202109*.nc')
#o.add_reader([reader_current, reader_wind, reader_wave]) 
o.add_reader([reader_current])
#o.set_config('environment:constant', {'x_sea_water_velocity': 0, 'y_sea_water_velocity': 1})
o.set_config('environment:constant', {'x_wind': 0, 'y_wind': 0})
##AO o.set_config("drift:max_age_seconds", 86401)

### Run Normal
#o.seed_elements(lon= -57,lat= 69,time=datetime(2021,9,1,12,00),number=100,radius=500,sail=10,draft=50,length=90,width=40)
#o.run(duration=timedelta(days=2))
##o.run(duration=timedelta(days=5),time_step=900, time_step_output=3600) #, outfile='test.nc')
#o.plot(fast=True)
##o.plot_property('draft')
'''

### Ensemble Testing : Randomise Ca and Co
# Define the number of ensemble members
nb_ensemble_members = 100
nb_particles=100

Ca_min=0.5 ; Ca_max=2.5
Co_min=0.5 ; Co_max=2.5

lon = -57.391907
lat = 72.082367
start_time = datetime(2021, 9, 1, 12, 00)


def generate_random_coeff(min_val, max_val, size, distribution_type='uniform'):
     if distribution_type == 'uniform':
          return np.random.uniform(min_val, max_val, size)
     elif distribution_type == 'normal':
          mean = (min_val + max_val) / 2
          std = (max_val - min_val) / 6  # 99.7% should fall within the range
          return np.clip(np.random.normal(mean, std, size), min_val, max_val)
     else:
          raise ValueError("Invalid distribution_type. Choose 'uniform' or 'normal'.")
          

### Run Ensemble
for i in range(nb_ensemble_members):
    # Define Output directory
    output_dir = 'Ens_results'
    os.makedirs(output_dir, exist_ok=True)

    o = OpenBerg()
    # Set configuration parameters
    o.set_config('drift:wave_rad', False)
    o.set_config('drift:stokes_drift', False)
    o.set_config('drift:coriolis', False)
    o.set_config('drift:vertical_profile', True)
    o.set_config('processes:grounding', False)
    o.set_config('processes:roll_over', False)
    # o.set_config('processes:melting', False)

    # Add reader for currents (and possibly winds, waves)
    reader_current = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/TP5Metno_m09_m10_y2021/202109*_dm-metno-MODEL-topaz5-ARC-b202109*-fv02.0.nc')
    reader_wind = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/Wind/6h.10m_wind_2021.nc')
    # reader_wave = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/Waves/MyWam3km_hindcast-cmems_202109*.nc')
    # o.add_reader([reader_current, reader_wind, reader_wave])
    o.add_reader([reader_current])

    # Set a constant wind environment
    #o.set_config('environment:constant', {'x_wind': 0, 'y_wind': 0})

    Co_rand_base = generate_random_coeff(Co_min, Co_max, nb_ensemble_members, distribution_type='normal')
    Ca_rand_base = generate_random_coeff(Ca_min, Ca_max, nb_ensemble_members, distribution_type='normal')
    
    # Generate random drag coefficients for each ensemble member
    Co_rand = Co_rand_base[i]
    Ca_rand = Ca_rand_base[i]
    #print("Ca_rand:", Ca_rand)
    #print("Co_rand:", Co_rand)
#    Ca_rand = Ca_rand_base.repeat(nb_particles / nb_ensemble_members)
#    Co_rand = Co_rand_base.repeat(nb_particles / nb_ensemble_members)
    print(f"Ensemble Member {i+1}: Water Drag Coefficient: {Co_rand}, Air Drag Coefficient: {Ca_rand}")

    # Seed elements with random drag coefficients
    o.seed_elements(lon=lon, lat=lat, time=start_time, number=nb_particles, radius=500,
                    water_drag_coeff=Co_rand, wind_drag_coeff=Ca_rand)

    # Run
    output_filename=os.path.join(output_dir, f'Ens_member_{i+1}.nc')
    o.run(duration=timedelta(days=3), time_step=900, time_step_output=3600, outfile=output_filename)
    #o.plot(linecolor='water_drag_coeff', buffer=.1)


# # Plots
simulations = []
for i in range(nb_ensemble_members):
     output_dir = 'Ens_results'
     os.makedirs(output_dir, exist_ok=True)
     output_filename=os.path.join(output_dir, f'Ens_member_{i+1}.nc')
     simulations.append(opendrift.open(output_filename))
colors = plt.cm.get_cmap('viridis', nb_ensemble_members)
# Assign each simulation a unique color
simulations[0].plot_comparison_colors = [colors(i) for i in range(nb_ensemble_members)]
# Create animation comparing the first simulation with the rest
simulations[0].animation(compare=simulations[1:], color='water_drag_coeff', filename=os.path.join(output_dir, 'Ens_animation.mp4'))







'''
#     o.seed_ensemble_v1(lon= -57,lat= 69, distribution_type='uniform', radius=100, number=300,
#                     time=start_time,
#                     sail=10,draft=50,length=90,width=40)
#     o.run(duration=timedelta(days=3),time_step=900, time_step_output=3600, outfile='Ens_Uniform.nc')
#     o.plot(fast=True)
#     o.animation(filename='Ens_Uniform_animation.mp4')

##############################################################
##############################################################
##############################################################
### Another Approach: All in one seed ###
Ca_min=0.5 ; Ca_max=2.5
Co_min=0.5 ; Co_max=2.5
ens_size = num_coeff =30 ; nb_particles_per_coeff=100      #(num_particles = num_coeff * nb_particles_per_coeff)
distribution_type='uniform'

def generate_random_coeff(min_val, max_val, size, distribution_type='uniform'):
        if distribution_type == 'uniform':
            return np.random.uniform(min_val, max_val, size)
        elif distribution_type == 'normal':
             mean = (min_val + max_val) / 2
             std = (max_val - min_val) / 6  # 99.7% should fall within the range
             return np.clip(np.random.normal(mean, std, size), min_val, max_val)
        else:
             raise ValueError("Invalid distribution_type. Choose 'uniform' or 'normal'.")


### Each number in the output array is independently drawn from a uniform distribution over the interval [min_val, max_val).
### This means that every number in this range has an equal chance of being chosen.

Ca_rand = generate_random_coeff(Ca_min, Ca_max, ens_size, distribution_type='normal').repeat(nb_particles_per_coeff)  
Co_rand = generate_random_coeff(Co_min, Co_max, ens_size, distribution_type='normal').repeat(nb_particles_per_coeff)


### The output array contains ens_size numbers that are evenly spaced between Ca_min and Ca_max.
### This means that the numbers are not random but are instead spaced at regular intervals.
##Ca_rand = np.linspace(Ca_min, Ca_max, ens_size)
##Co_rand = np.linspace(Co_min, Co_max, ens_size)

o.seed_elements(lon=lon, lat=lat, radius=1000, time=start_time, number=len(Co_rand),
                wind_drag_coeff=Ca_rand, water_drag_coeff=Co_rand)

o.run(duration=timedelta(days=20),time_step=900, time_step_output=3600, outfile='Ens_Run.nc')
o.plot(fast=True)
o.animation(filename='Ens_animation.mp4', color=wind_drag_coeff))


o.animationr(color=wind_drag_coeff)
##############################################################
##############################################################
##############################################################
'''
####################################################################################################################################

###### Iceberg 1 West Greenland 2021 ### 11_DSE_300434065758620
#lon = -56.391907
#lat = 72.082367
#start_time = datetime(2021,9,22,14,32)
#o.seed_elements(lon,lat,radius=500,number=100,time=start_time,sail=45,draft=220,length=220,width=420)

###### Iceberg 2 Fram Strait 2018
#lon = -20.7183647984
#lat = 72.5984112684
#start_time = datetime(2018,7,10,18,57)
#o.seed_elements(lon,lat,radius=500,number=100,time=start_time,sail=35,draft=175,length=1500,width=700)

###### Iceberg 3 Iceberg  West Greenland 2019 
#lon = -76.0735
#lat = 77.18793
#start_time = datetime(2019,8,1,12,17,29)
#o.seed_elements(lon,lat,radius=500,number=100,time=start_time,sail=20,draft=100,length=2600,width=2400)

###### Iceberg 4 Iceberg  West Greenland 2020 
#lon = -77.6438
#lat = 77.08413
#start_time = datetime(2020,8,1,12,17,19) 
#o.seed_elements(lon,lat,radius=500,number=100,time=start_time,sail=20,draft=100,length=2000,width=1700)


''' Test 

# Define colors for each ensemble member (customize as needed)
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']  # Add more colors if needed

# Create a figure for the combined plot (optional, for visualization)
plt.figure()

# Create a list to store each ensemble member's object for animation
ensemble_members_objects = []

for i in range(nb_ensemble_members):
    o = OpenBerg()
    o.set_config('drift:wave_rad', False)
    o.set_config('drift:stokes_drift', False)
    o.set_config('drift:coriolis', False)
    o.set_config('drift:vertical_profile', True)
    o.set_config('processes:grounding', False)
    o.set_config('processes:roll_over', False)
    # o.set_config('processes:melting', False)
    
    reader_current = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/TP5Metno_m09_m10_y2021/202109*_dm-metno-MODEL-topaz5-ARC-b2021*-fv02.0.nc') ###TP5 Metno 3D
    #reader_wind = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/Wind/6h.10m_wind_2021.nc')
    #reader_wave = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/Waves/MyWam3km_hindcast-cmems_202109*.nc')
    #o.add_reader([reader_current, reader_wind, reader_wave]) 
    o.add_reader([reader_current])
    o.set_config('environment:constant', {'x_wind': 0, 'y_wind': 0})

    ##o.add_reader([reader_current])
    ##o.set_config('environment:constant', {'x_wind': 0, 'y_wind': 0})

    # Generate random drag coefficients for each ensemble member
    Ca_rand = generate_random_coeff(Ca_min, Ca_max, nb_ensemble_members, distribution_type='normal').repeat(nb_particles/nb_ensemble_members)  
    Co_rand = generate_random_coeff(Co_min, Co_max, nb_ensemble_members, distribution_type='normal').repeat(nb_particles/nb_ensemble_members)  

    print(f"Ensemble Member {i+1}: Water Drag Coefficient: {Co_rand}, Air Drag Coefficient: {Ca_rand}")

    # Seed elements with random drag coefficients
    o.seed_elements(lon=lon, lat=lat, time=start_time, number=nb_particles, radius=500,
                    water_drag_coeff=Co_rand, wind_drag_coeff=Ca_rand)

    # Run simulation
    o.run(duration=timedelta(days=3), time_step=900, time_step_output=3600, outfile='Ens_Test.nc')
    


    # Generate and save the plot for each ensemble member with a unique color
    color = colors[i % len(colors)]  # Assign a color for the current ensemble member
    result_file = os.path.join(Output_dir, f'Fig_ensemble_member_{i+1}.png')
    o.plot(fast=True, filename=result_file, color=lon) #### AO

    ##Save the draft property plot for the ensemble member
    ##draft_file = os.path.join(Output_dir, f'Fig_draft_ensemble_member_{i+1}.png')
    ##o.plot_property('draft', filename=draft_file, color=color)

    # Store the ensemble member object for creating a combined animation
    ensemble_members_objects.append(o)

# After the loop, create a combined animation with all ensemble members in different colors
animation_file = os.path.join(Output_dir, 'Ens_Normal_animation_all_members.mp4')
for o in ensemble_members_objects:
    o.animation(filename=animation_file)  # Animate all the stored ensemble members
'''
####################################################################################################################################


