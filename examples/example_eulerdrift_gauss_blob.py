"""
Euler simulation / Finite differnce of Gaussian blob
=====================================================
"""
import matplotlib.pyplot as plt
from opendrift.models import eulerdrift

eulerdrift.install_logs() #remove

s = eulerdrift.ExplSimulation.new(shape = (200, 200))
s.readers.append(
    eulerdrift.readers.ConstantReader.new_xy(5, 5))

#%%
# Set the diffusion, greater diffusion increases stability of simulation.
s.D = 10.

#%%
# Add a source field (Gaussian blob) at the center of the grid
loc, lac = s.grid.center()
s.source_gaussian_blob(loc, lac, 1., 70, 20.)

#%%
# Plot the initial conditions
s.grid.plot()
plt.show()

#%%
# Integrate
s.integrate(dt = .01, max_steps = 5000)

#%%
# Plot integrated field
s.grid.plot()
plt.show()
