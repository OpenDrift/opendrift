#!/usr/bin/env python
"""
ChemicalDrift - Transport and fate of organic compounds
========================================================
"""

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.chemicaldrift import ChemicalDrift
from opendrift.readers.reader_constant import Reader as ConstantReader
from datetime import timedelta, datetime
import matplotlib.pyplot as plt
import numpy as np


o = ChemicalDrift(loglevel=0, seed=0)

# Norkyst
reader_norkyst = reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
mixed_layer = ConstantReader({'ocean_mixed_layer_thickness': 40})
o.add_reader([reader_norkyst,mixed_layer])

o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Large1994')
o.set_config('chemical:particle_diameter',30.e-6)  # m
o.set_config('chemical:particle_diameter_uncertainty',5.e-6) # m
o.set_config('chemical:sediment:resuspension_critvel',0.15) # m/s
o.set_config('chemical:transformations:volatilization', True)
o.set_config('chemical:transformations:degradation', True)
o.set_config('chemical:transformations:degradation_mode', 'OverallRateConstants')
o.set_config('seed:LMM_fraction',.9)
o.set_config('seed:particle_fraction',.1)
o.set_config('general:coastline_action', 'previous')

o.init_chemical_compound("Phenanthrene")

# Modify half-life times with unrealistic values for this demo
o.set_config('chemical:transformations:t12_W_tot', 6.) # hours
o.set_config('chemical:transformations:t12_S_tot', 12.) # hours

o.list_configspec()


#%%
# Seeding 500 lagrangian elements each representign 2mg og target chemical

td=datetime.today()
time = td - timedelta(days=10)

latseed= 57.6;   lonseed= 10.6 # Skagen

ntraj=500
iniz=np.random.rand(ntraj) * -10. # seeding the chemicals in the upper 10m

o.seed_elements(lonseed, latseed, z=iniz, radius=2000,number=ntraj,time=time, mass=2e3) 


#%%
# Running model

o.run(steps=48*2, time_step=1800, time_step_output=1800)


#%%
# Print and plot results

print('Final speciation:')
for isp,sp in enumerate(o.name_species):
    print ('{:32}: {:>6}'.format(sp,sum(o.elements.specie==isp)))

print('Number of transformations:')
for isp in range(o.nspecies):
    print('{}'.format(['{:>9}'.format(np.int32(item)) for item in o.ntransformations[isp,:]]) )

m_pre = o.result.mass.isel(time=-1).sum()
m_deg = o.result.mass_degraded.isel(time=-1).sum()
m_vol = o.result.mass_volatilized.isel(time=-1).sum()
m_tot = m_pre + m_deg + m_vol

print('Mass budget for target chemical:')
print('mass preserved       : {:.3f}'.format(m_pre * 1e-6),' g   {:.3f}'.format(m_pre/m_tot*100),'%')
print('mass degraded        : {:.3f}'.format(m_deg * 1e-6),' g   {:.3f}'.format(m_deg/m_tot*100),'%')
print('mass volatilized     : {:.3f}'.format(m_vol * 1e-6),' g   {:.3f}'.format(m_vol/m_tot*100),'%')

legend=['dissolved', '', 'SPM', 'sediment', '']

o.animation_profile(color='specie',
            markersize='mass',
            markersize_scaling=30,
            alpha=.5,
            vmin=0,vmax=o.nspecies-1,
            legend = legend,
            legend_loc = 3,
            fps = 10
            )
#%%
# .. image:: /gallery/animations/example_chemicaldrift_0.gif

o.animation(color='specie',
            markersize='mass',
            markersize_scaling=100,
            alpha=.5,
            vmin=0,vmax=o.nspecies-1,
            colorbar=False,
            fast = True,
            legend = legend,
            legend_loc = 3,
            fps=10
            )
#%%
# .. image:: /gallery/animations/example_chemicaldrift_1.gif


#%%

o.plot_mass(legend = legend,
            time_unit = 'hours',
            title = 'Chemical mass budget')
