"""
Radionuclides
=============
"""
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
reader_norkyst = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')

#o.add_reader([reader_norkyst, reader_arome])
o.add_reader([reader_norkyst])

#time = reader_norkyst.start_time
td=datetime.today()
time = datetime(td.year, td.month, td.day,0)

#latseed= 61.2; lonseed= 4.3    # Sognesjen
latseed= 59.0;   lonseed= 10.75 # Hvaler/Koster



ntraj=2000
iniz=np.random.rand(ntraj) * -10.         # seeding the radionuclides in the upper 10m
init_speciation = np.ones(ntraj)*0
diam=np.zeros_like(init_speciation,dtype=np.float32)
o.seed_elements(lonseed, latseed, z=iniz, radius=100,number=ntraj,
                time=time,
                diameter=diam, specie=init_speciation)

init_speciation = np.ones(ntraj)*1
diam=np.zeros_like(init_speciation,dtype=np.float32) * 5.e-6
o.seed_elements(lonseed, latseed, z=iniz, radius=100,number=ntraj,
                time=time,
                diameter=diam, specie=init_speciation)

# init_speciation = np.ones(ntraj)*2
# o.seed_elements(lonseed, latseed, z=iniz, radius=100,number=ntraj,
#                 time=time,
#                 diameter=diam, specie=init_speciation)

# Adjusting some configuration
o.set_config('processes:turbulentmixing', True)
o.set_config('turbulentmixing:diffusivitymodel','zero')  # include settling without vertical turbulent mixing
#o.set_config('turbulentmixing:diffusivitymodel','environment')  # include settling without vertical turbulent mixing
# Vertical mixing requires fast time step
o.set_config('turbulentmixing:timestep', 600.) # seconds

# Activate the desired species
o.set_config('radionuclide:species:LMM', True)
#o.set_config('radionuclide:species:LMMcation', True)
#o.set_config('radionuclide:species:LMManion', True)
#o.set_config('radionuclide:species:Colloid', True)
#o.set_config('radionuclide:species:Humic_colloid', True)
#o.set_config('radionuclide:species:Polymer', True)
o.set_config('radionuclide:species:Particle_reversible', True)
#o.set_config('radionuclide:species:Particle_irreversible', True)
o.set_config('radionuclide:species:Particle_slowly_reversible', True)

o.set_config('radionuclide:species:Sediment_reversible', True)
#o.set_config('radionuclide:species:Sediment_irreversible', True)
o.set_config('radionuclide:species:Sediment_slowly_reversible', True)

o.set_config('radionuclide:particle_diameter',5.e-6)  # m
#o.set_config('radionuclide:particle_diameter_uncertainty',1.e-7) # m
o.set_config('radionuclide:transformations:Kd',2.e0) # (m3/kg)
#o.set_config('radionuclide:transformations:slow_coeff',1.e-6)

o.set_config('radionuclide:sediment:resuspension_depth',2.)
o.set_config('radionuclide:sediment:resuspension_depth_uncert',0.1)
o.set_config('radionuclide:sediment:resuspension_critvel',0.15)


#
#o.set_config('radionuclide:transfer_setup','dummy')
o.set_config('radionuclide:transfer_setup','Bokna_137Cs')
#o.set_config('radionuclide:transfer_setup','Sandnesfj_Al')

# By default, radionuclides do not strand towards coastline
o.set_config('general:coastline_action', 'previous')

#o.set_config('general:use_auto_landmask',False)

o.list_configspec()

# Running model (until end of driver data)
o.run(steps=48*2, time_step=1800,outfile='radio.nc')


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
            filename='radionuclides.gif'

#            fast = True
            )
#%%
# .. image:: /gallery/animations/radionuclides.gif

#o.plot_vertical_distribution()
#o.plot_property('specie')
o.animation_profile(filename='radionuclides_profile.gif')
#%%
# .. image:: /gallery/animations/radionuclides_profile.gif

o.plot(linecolor='specie',vmin=0,vmax=o.nspecies-1,#fast=True,
#       background='land_binary_mask',
       #background='sea_floor_depth_below_sea_level'
       )


