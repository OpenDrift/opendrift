import opendrift
import matplotlib.pyplot as plt 
import os

nb_ensemble_members = 100


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
simulations[0].animation(compare=simulations[1:], markersize='water_drag_coeff', buffer=.01, filename=os.path.join(output_dir, 'Ens_animation.mp4')) #,linewidth=20

#o.plot(fast=True, compare=o2, legend=['Current + 3 % wind drift', 'Current only'])
#o.animation(fast=True, compare=o2, legend=['Current + 3 % wind drift', 'Current only'])