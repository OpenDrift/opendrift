
#!/usr/bin/env python
"""
Icebergs (openberg)
====================
************************ Data sources *********************
*. TOPAZ5: 2021-07-05  to  2024-03-21        cmems_mod_arc_phy_anfc_6km_detided_P1D-m     ###or cmems_mod_arc_phy_anfc_6km_detided_PT1H-i
*. TOPAZ6: 2018-01-01  to  2024-03-21        dataset-topaz6-arc-15min-3km-be              ###2D tidal model
*. TOPAZ4: 1991-01-01  to  2022-12-31        cmems_mod_arc_phy_my_topaz4_P1D-m            ###No more used
*. GLORYS: 1993-01-01  to  2022-07-31        cmems_mod_glo_phy_myint_0.083deg_P1D-m       ###or  cmems_mod_glo_phy_my_0.083deg_P1D-m
*. Globcurrent: 1993-01-01  to  2022-07-31   cmems_obs_mob_glo_phy-cur_my_0.25deg_P1D-m   ###Sat product
*. Waves: 1977-01-01  to  2022-11-30         cmems_mod_arc_wav_my_3km_PT1H-i
"""

#import copernicus_marine_client as copernicusmarine
from opendrift.readers.reader_netCDF_CF_generic import Reader
from opendrift.models.openberg import OpenBerg
from datetime import datetime,timedelta
import numpy as np
import os

### Load credentials of CMEMS Account from an external file
#with open("CMEMS_credentials.txt") as f:
#     lines = f.read().splitlines()
#     if len(lines) >= 2:
#          username = lines[0]
#          password = lines[1]

#USERNAME=username
#PASSWORD=password

### Define datasets
#Wind = '/Volumes/ACH_HDD/DATA_ACCIBERG/Wind/6h.10m_wind_2021.nc'
#Ocean = 'cmems_mod_glo_phy_myint_0.083deg_P1D-m'
#Wave = 'cmems_mod_glo_wav_my_0.2deg_PT3H-i'

### Add Readers
#reader_wind = Reader(Wind)

#ds_current = copernicusmarine.open_dataset(dataset_id=Ocean, username=USERNAME, password=PASSWORD)
#reader_current = Reader(ds_current)

#ds_wave = copernicusmarine.open_dataset(dataset_id=Wave, username=USERNAME, password=PASSWORD)
#reader_wave = Reader(ds_wave)

###################################################################
### Run Normal
o = OpenBerg()
# Set configuration parameters
o.set_config('drift:stokes_drift', False)
o.set_config('drift:vertical_profile', True)
o.set_config('drift:wave_rad', False)
o.set_config('drift:coriolis', False)
o.set_config('drift:sea_surface_slope', False)
o.set_config('processes:grounding', False)
o.set_config('processes:roll_over', False)
o.set_config('processes:melting', False)

reader_current = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/TP5Metno_m09_m10_y2021/202109*_dm-metno-MODEL-topaz5-ARC-b2021*-fv02.0.nc') ###TP5 Metno 3D
#reader_wind = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/Wind/6h.10m_wind_2021.nc')
#reader_wave = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/Waves/MyWam3km_hindcast-cmems_202109*.nc')
o.add_reader([reader_current])

o.seed_elements(lon= -57,lat= 69,time=datetime(2021,9,1,12,00),number=1,radius=500,sail=10,draft=50,length=90,width=40)
o.run(duration=timedelta(days=7),time_step=900, time_step_output=3600, outfile='Test.nc')
o.plot(fast=True)

# ###################################################################
# ### Run Ensembles
# ### Define parameters needed for the run: Ensemble Testing with random Ca and Co coeff.
# nb_ensemble_members = 30  # number of ensemble members
# nb_particles = 100         # number of particles for each member
# Ca_min = 0.5 ; Ca_max = 2.5
# Co_min = 0.5 ; Co_max = 2.5
# lon = -57.391907 # starting point for the seeding
# lat = 72.082367
# start_time = datetime(2021, 9, 1, 12, 00)

# ### Function to generate random coefficients from unifrom or normal distribution
# def generate_random_coeff(min_val, max_val, size, distribution_type='uniform'):
#      if distribution_type == 'uniform':
#           return np.random.uniform(min_val, max_val, size)
#      elif distribution_type == 'normal':
#           mean = (min_val + max_val) / 2
#           std = (max_val - min_val) / 6  # 99.7% should fall within the range
#           return np.clip(np.random.normal(mean, std, size), min_val, max_val)
#      else:
#           raise ValueError("Invalid distribution_type. Choose 'uniform' or 'normal'.")
          
# ### Run Ensemble
# for i in range(nb_ensemble_members):
#     # Define Output directory
#     output_dir = 'Ens_results'
#     os.makedirs(output_dir, exist_ok=True)

#     o = OpenBerg()
#     # Set configuration parameters
#     o.set_config('drift:wave_rad', False)
#     o.set_config('drift:stokes_drift', False)
#     o.set_config('drift:coriolis', False)
#     o.set_config('drift:sea_surface_slope', False) #### Test
#     o.set_config('drift:vertical_profile', True)
#     o.set_config('processes:grounding', False)
#     o.set_config('processes:roll_over', False)
#     # o.set_config('processes:melting', False)

#     # Add reader for currents (possibly winds, waves)
#     reader_current = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/TP5Metno_m09_m10_y2021/202109*_dm-metno-MODEL-topaz5-ARC-b202109*-fv02.0.nc')  # from local machine
#     reader_wind = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/Wind/6h.10m_wind_2021.nc')  # from local machine
#     #reader_wave = Reader('/Volumes/ACH_HDD/DATA_ACCIBERG/Waves/MyWam3km_hindcast-cmems_202109*.nc')  # from local machine
#     # o.add_reader([reader_current, reader_wind, reader_wave])
#     o.add_reader([reader_current, reader_wind])

#     # Set a constant wind environment
#     #o.set_config('environment:constant', {'x_wind': 0, 'y_wind': 0})

#     # Generate random drag coefficients for each ensemble member
#     Co_rand_base = generate_random_coeff(Co_min, Co_max, nb_ensemble_members, distribution_type='uniform')
#     Ca_rand_base = generate_random_coeff(Ca_min, Ca_max, nb_ensemble_members, distribution_type='uniform')
#     Co_rand = Co_rand_base[i]
#     Ca_rand = Ca_rand_base[i]
#     print(f"Ensemble Member {i+1}: Water Drag Coefficient: {Co_rand}, Air Drag Coefficient: {Ca_rand}")

#     # Seed elements with random drag coefficients
#     # Latitude and longitude define the starting point for the seeding (initialisation) of particles.
#     # The radius is an uncertainty radius (standard deviation) around the given location within which all the elements will be seeded.
#     o.seed_elements(lon=lon, lat=lat, time=start_time, number=nb_particles, radius=500,
#                     water_drag_coeff=Co_rand, wind_drag_coeff=Ca_rand)

#     # Run
#     output_filename=os.path.join(output_dir, f'Ens_member_{i+1}.nc')
#     o.run(duration=timedelta(days=21), time_step=900, time_step_output=3600, outfile=output_filename)
# ###################################################################