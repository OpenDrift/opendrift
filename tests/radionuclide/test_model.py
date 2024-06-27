from . import *
from opendrift.readers import reader_netCDF_CF_generic, reader_ROMS_native
from opendrift.models.radionuclides import RadionuclideDrift
from datetime import timedelta, datetime
import numpy as np

def test_radio_nuclide(test_data):
    o = RadionuclideDrift(loglevel=20, seed=0)  # Set loglevel to 0 for debug information

    # Norkyst
    reader_norkyst = reader_netCDF_CF_generic.Reader(test_data + '/14Jan2016_NorKyst_z_3d/NorKyst-800m_ZDEPTHS_his_00_3Dsubset.nc')
    o.add_reader([reader_norkyst])

    # Adjusting some configuration
    o.set_config('drift:vertical_mixing', True)
    o.set_config('vertical_mixing:diffusivitymodel','environment')  # include settling without vertical turbulent mixing
    # Vertical mixing requires fast time step
    o.set_config('vertical_mixing:timestep', 600.) # seconds
    o.set_config('drift:horizontal_diffusivity', 10)
    o.set_config('radionuclide:particle_diameter',5.e-6)  # m
    o.set_config('radionuclide:isotope','137Cs')
    o.set_config('radionuclide:specie_setup','LMM + Rev')

    o.set_config('general:coastline_action', 'previous')
    o.set_config('seed:LMM_fraction',.45)
    o.set_config('seed:particle_fraction',.55)

    time = reader_norkyst.start_time

    latseed=62.0;
    lonseed=4.3 # Bergen (?)

    ntraj=5000
    iniz=np.random.rand(ntraj) * -10. # seeding the radionuclides in the upper 10m

    o.seed_elements(lonseed, latseed, z=iniz, radius=1000,number=ntraj,
                    time=time,)


    o.run(steps=4, time_step=1800, time_step_output=3600)
