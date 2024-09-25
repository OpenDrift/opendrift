import opendrift
import matplotlib.pyplot as plt 
import os

nb_ensemble_members = 30
plt.rcParams['figure.figsize'] = (30, 20)
#plt.rcParams['legend.fontsize'] = 'small'


# Plots
simulations = []
for i in range(nb_ensemble_members):
    output_dir = 'Ens_results'
    os.makedirs(output_dir, exist_ok=True)
    output_filename=os.path.join(output_dir, f'Ens_member_{i+1}.nc')
    simulations.append(opendrift.open(output_filename))

colors= plt.cm.get_cmap('jet', nb_ensemble_members)
simulations[0].plot_comparison_colors = [colors(i) for i in range(nb_ensemble_members)]
#legend_labels = [f'Member_{i+1}' for i in range(nb_ensemble_members)]

# Create animation comparing the first simulation with the rest
output_animation_filename = os.path.join(output_dir, 'Ens_animation.mp4')
simulations[0].animation(compare=simulations[1:], fast=False, buffer=.01, filename=output_animation_filename) # legend=legend_labels,
