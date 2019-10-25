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
#reader_bokna = reader_ROMS_native.Reader('/home/magnes/projects/MARTINI/data/norfjords_160m_his.nc4_2015110?01-2015110?00')

#o.add_reader([reader_norkyst, reader_arome])
o.add_reader([reader_norkyst])
#o.add_reader([reader_bokna])



# release radionuclides at defined position and time
#time = reader_norkyst.start_time
td=datetime.today()
time = datetime(td.year, td.month, td.day,0)
#time = datetime(2015,11,2,1)   # Boknafj

latseed= 61.2; lonseed= 4.3    # Sognesjen
#latseed= 59.0;   lonseed= 10.75 # Hvaler/Koster
#latseed= 59.486; lonseed= 5.883   # Boknafjorden, Vikedal
#latseed= 59.09; lonseed= 5.39   # Boknafjorden, Kvitsoy


ntraj=10000
#init_speciation = np.random.randint(o.nspecies,size=ntraj)
init_speciation = np.ones(ntraj)*0

diam=np.zeros_like(init_speciation,dtype=np.float32)
#diam[init_speciation==o.name_species.index('Particle reversible')] = 0.001
#diam[init_speciation==o.name_species.index('Particle irreversible')] = 0.001
#print 'diam',diam


o.seed_elements(lonseed, latseed, z=-4., radius=1000,number=ntraj,
                time=time, diameter=diam, specie=init_speciation)



# Adjusting some configuration
o.set_config('processes:turbulentmixing', True)
o.set_config('turbulentmixing:diffusivitymodel','zero')  # include settling without vertical turbulent mixing
# Vertical mixing requires fast time step
o.set_config('turbulentmixing:timestep', 600.) # seconds

# Activate the desired species
o.set_config('radionuclide:species:LMM', True)
#o.set_config('radionuclide:species:Colloid', True)
o.set_config('radionuclide:species:Particle_reversible', True)
#o.set_config('radionuclide:species:Particle_irreversible', True)
o.set_config('radionuclide:species:Particle_slowly_reversible', True)

o.set_config('radionuclide:species:Sediment_reversible', True)
#o.set_config('radionuclide:species:Sediment_irreversible', True)
o.set_config('radionuclide:species:Sediment_slowly_reversible', True)

#o.set_config('radionuclide:particle_diameter',100.e-6)  # m
#o.set_config('radionuclide:particle_diameter_uncertainty',1.e-7) # m
o.set_config('radionuclide:transformations:Kd',2.e1) # (m3/kg)
o.set_config('radionuclide:transformations:slow_coeff',1.e-6)

#
#o.set_config('radionuclide:transfer_setup','dummy')
o.set_config('radionuclide:transfer_setup','Bokna_137Cs')

#o.set_config('general:use_basemap_landmask',False)

o.list_configspec()

# Running model (until end of driver data)
o.run(steps=4*2, time_step=1800,outfile='radio.nc')

# Print and plot results
print(o)
print('Final speciation:')
for isp,sp in enumerate(o.name_species):
    print ('{:32}: {:>6}'.format(sp,sum(o.elements.specie==isp)))


o.animation(color='specie',
            vmin=0,vmax=o.nspecies,
            colorbar=True)
#o.plot_vertical_distribution()
#o.plot_property('specie')
o.animation_profile()
#o.animation()
o.plot(linecolor='specie')
#o.animation(density=True,unitfactor=0.1,vmax=1400)

