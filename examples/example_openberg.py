#!/usr/bin/env python
"""
Icebergs (openberg)
====================
"""

from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import opendrift
from opendrift.models.openberg import OpenBerg

#%%
# Currents and wind forcing
forcing = ['https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be',
           'https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd']

#%%
# A permutation of iceberg sizes/dimensions
n = 10
lengths = np.linspace(50, 500, n)
widths = np.linspace(20, 200, n)
sails = np.linspace(5, 50, n)
drafts = np.linspace(2, 120, n)
lengths, widths, sails, drafts = np.meshgrid(lengths, widths, sails, drafts)

icebergs = {'lon': 18.127, 'lat': 74.776, 'time': datetime.now(),
            'number': lengths.size, 'radius': 500,
            'sail': sails, 'draft': drafts, 'length': lengths, 'width': widths}

#%%
# Simulating drift for 48 hours
o = OpenBerg()
o.set_config('drift:vertical_profile', False)
o.set_config('drift:horizontal_diffusivity', 100)
o.add_readers_from_list(forcing)
o.seed_elements(**icebergs)
o.run(duration=timedelta(hours=48))

o.animation(color='draft', contourlines=np.arange(0, 500, 25))

#%%
# .. image:: /gallery/animations/example_openberg_0.gif

o.plot(contourlines=np.arange(0, 500, 25))

#%%
# Plotting the speed of icebergs
iceberg_speed = np.sqrt(o.result.iceb_x_velocity**2 + o.result.iceb_y_velocity**2)
iceberg_speed.plot.line(x='time', add_legend=False, color='gray')
plt.ylabel('Iceberg speed  [m/s]')
plt.show()