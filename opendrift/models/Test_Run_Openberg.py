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
    o.add_reader([reader_current, reader_wind])

    # Set a constant wind environment
    #o.set_config('environment:constant', {'x_wind': 0, 'y_wind': 0})

    Co_rand_base = generate_random_coeff(Co_min, Co_max, nb_ensemble_members, distribution_type='normal')
    Ca_rand_base = generate_random_coeff(Ca_min, Ca_max, nb_ensemble_members, distribution_type='normal')
    
    # Generate random drag coefficients for each ensemble member
    Co_rand = Co_rand_base[i]
    Ca_rand = Ca_rand_base[i]
    print(f"Ensemble Member {i+1}: Water Drag Coefficient: {Co_rand}, Air Drag Coefficient: {Ca_rand}")

    # Seed elements with random drag coefficients
    o.seed_elements(lon=lon, lat=lat, time=start_time, number=nb_particles, radius=500,
                    water_drag_coeff=Co_rand, wind_drag_coeff=Ca_rand)

    # Run
    output_filename=os.path.join(output_dir, f'Ens_member_{i+1}.nc')
    o.run(duration=timedelta(days=3), time_step=900, time_step_output=3600, outfile=output_filename)


# Plots
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
