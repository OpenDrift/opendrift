#!/usr/bin/env python


from opendrift.readers import reader_netCDF_CF_generic, reader_ROMS_native
from opendrift.models.radionuclides import RadionuclideDrift
from datetime import timedelta, datetime
import numpy as np
from numpy import dtype




o = RadionuclideDrift(
                      loglevel=20,seed=0)  # Set loglevel to 0 for debug information

# Arome
#reader_arome = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')
#reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/meps25files/meps_det_extracted_2_5km_latest.nc')

# Norkyst
#reader_norkyst = reader_netCDF_CF_generic.Reader(o.test_data_folder() + '/14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc') 
#reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

#o.add_reader([reader_norkyst, reader_arome])
#o.add_reader([reader_norkyst])

#time = reader_norkyst.start_time
td=datetime.today()
#time = datetime(td.year, td.month, td.day,0)

#latseed= 61.2; lonseed= 4.3    # Sognesjen
#latseed= 59.0;   lonseed= 10.75 # Hvaler/Koster

#  Boknafjorden case
#reader_bokna = reader_ROMS_native.Reader('/home/magnes/projects/openDrift-RAD/data/norfjords_160m_his.nc4_2015110?01-2015110?00')
#o.add_reader([reader_bokna])
#time = datetime(2019,9,2,1)   # Boknafj
#time = datetime(2015,11,2,1)   # Boknafj
#latseed= 59.486; lonseed= 5.883   # Boknafjorden, Vikedal
#latseed= 59.09; lonseed= 5.39   # Boknafjorden, Kvitsoy

# Sandnesfjorden case
#reader_global = reader_global_landmask.Reader()
#reader_sandnesfj = reader_ROMS_native.Reader('/home/magnes/projects/openDrift-RAD/data/norfjords_32m_his.nc4_20080???01-2008????00')
#o.add_reader([reader_sandnesfj])
#time = datetime(2008,4,21,1)   # Sandnesfj
#latseed= 58.6872; lonseed= 9.0919  # Laget, Sandnesfjorden, Risr
#latseed= 58.6825; lonseed= 9.0717  # Laget, Sandnesfjorden, Risr


# Sandnesfjorden 2 case
reader_sandnesfj = reader_ROMS_native.Reader('/home/magnes/projects/openDrift-RAD/data/norfjords_32m_his.nc4_20180???00-2018????23')
o.add_reader([reader_sandnesfj])
time = datetime(2018,6,20,5)   # Sandnesfj
latseed= 58.6722; lonseed= 8.9870  # Storelva, Sandnesfjorden, Risr
#latseed= 58.6731; lonseed= 9.00479  # Sundet mellom Songev og Nevestadfj, Sandnesfj, Risr

#print reader_sandnesfj


ntraj=2000
#init_speciation = np.random.randint(o.nspecies,size=ntraj)
iniz=np.random.rand(ntraj) * -1.5
init_speciation = np.ones(ntraj)*0
diam=np.zeros_like(init_speciation,dtype=np.float32)
o.seed_elements(lonseed, latseed, z=iniz, radius=21,number=ntraj,
                #time=[time, time+timedelta(hours=4)],
                time=time,
                diameter=diam, specie=init_speciation)

init_speciation = np.ones(ntraj)*2
o.seed_elements(lonseed, latseed, z=iniz, radius=21,number=ntraj,
                time=time,
                #time=[time, time+timedelta(hours=24)],
                diameter=diam, specie=init_speciation)

init_speciation = np.ones(ntraj)*4
o.seed_elements(lonseed, latseed, z=iniz, radius=21,number=ntraj,
                time=time,
                #time=[time, time+timedelta(hours=24)],
                diameter=diam, specie=init_speciation)

# Adjusting some configuration
o.set_config('processes:turbulentmixing', True)
o.set_config('turbulentmixing:diffusivitymodel','zero')  # include settling without vertical turbulent mixing
#o.set_config('turbulentmixing:diffusivitymodel','environment')  # include settling without vertical turbulent mixing
# Vertical mixing requires fast time step
o.set_config('turbulentmixing:timestep', 200.) # seconds

# Activate the desired species
#o.set_config('radionuclide:species:LMM', True)
o.set_config('radionuclide:species:LMMcation', True)
o.set_config('radionuclide:species:LMManion', True)
#o.set_config('radionuclide:species:Colloid', True)
o.set_config('radionuclide:species:Humic_colloid', True)
o.set_config('radionuclide:species:Polymer', True)
o.set_config('radionuclide:species:Particle_reversible', True)
#o.set_config('radionuclide:species:Particle_irreversible', True)
#o.set_config('radionuclide:species:Particle_slowly_reversible', True)

o.set_config('radionuclide:species:Sediment_reversible', True)
#o.set_config('radionuclide:species:Sediment_irreversible', True)
#o.set_config('radionuclide:species:Sediment_slowly_reversible', True)

o.set_config('radionuclide:particle_diameter',5.e-6)  # m
#o.set_config('radionuclide:particle_diameter_uncertainty',1.e-7) # m
o.set_config('radionuclide:transformations:Kd',2.e1) # (m3/kg)
o.set_config('radionuclide:transformations:slow_coeff',1.e-6)

o.set_config('radionuclide:sediment:resuspension_depth',2.)
o.set_config('radionuclide:sediment:resuspension_depth_uncert',0.1)
o.set_config('radionuclide:sediment:resuspension_critvel',0.005)


#
#o.set_config('radionuclide:transfer_setup','dummy')
#o.set_config('radionuclide:transfer_setup','Bokna_137Cs')
o.set_config('radionuclide:transfer_setup','Sandnesfj_Al')

# By default, radionuclides do not strand towards coastline
o.set_config('general:coastline_action', 'previous')

#o.set_config('general:auto_landmask_resolution','l')
o.set_config('general:use_auto_landmask',False)

o.list_configspec()

# Running model (until end of driver data)
o.run(steps=2*9, time_step=400,outfile='../MyTests/radio.nc')

# Print and plot results
print(o)
print('Final speciation:')
for isp,sp in enumerate(o.name_species):
    print ('{:32}: {:>6}'.format(sp,sum(o.elements.specie==isp)))

print('Number of transformations:')
print (o.ntransformations)


o.animation(color='specie',
            vmin=0,vmax=o.nspecies-1,
            colorbar=True,
            fast = True
            )
#o.plot_vertical_distribution()
#o.plot_property('specie')
#o.animation_profile()
o.animation(fast=True,
            #background='land_binary_mask',
            #background='sea_floor_depth_below_sea_level'
            )
o.plot(linecolor='specie',vmin=0,vmax=o.nspecies-1,fast=True,
       #background='land_binary_mask',
       #background='sea_floor_depth_below_sea_level'
       )
#o.animation(density=True,unitfactor=0.1,vmax=1400)

