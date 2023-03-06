# %%
import numpy as np
from simulation import run_simulation
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as colors
from viz import get_patch_sizes

def plot_slices_of_simulation_biomass(results_dict, title = "", plot_burn_masks=True): 
    '''Params: 
    times: a list of time steps 
    slices: an array with shape (2*save_times, _, _), where L is the side length of a 2-D world. 
        Each subarray contains the state of the forest pre- and post-fire at a time in save_times
    burn_mask (optional): array of binary masks, if included, burned areas will be displayed in orange
    Graphs the state of the simulation in a series of subplots, showing its progression over time 
    '''
    slices = results_dict["output_slices"]
    times = results_dict["output_times"]
    tree_cc = results_dict["params_dict"]["tree_carrying_capacity"]
    if plot_burn_masks: 
        num_rows = 1
        burn_masks = results_dict["burned_area_masks"]
    else: 
        num_rows = 1

    fig, ax = plt.subplots(nrows=num_rows, ncols=len(times), figsize=(20, 8))
    for i in range(len(times)): 
        ax[i].xaxis.set_visible(False) #hide axis ticks 
        ax[i].yaxis.set_visible(False)
        ax[i].xaxis.set_visible(False) 
        ax[i].yaxis.set_visible(False)

        # show the plot after growth season 
        biomass_upper_bound = 1e6
        norm = colors.Normalize(vmin=0, vmax=biomass_upper_bound)
        ax[i].imshow(slices[2*i], vmin=1, vmax=tree_cc*1.2, cmap = "summer_r", alpha=(slices[2*i] > 0).astype(float))
        

        # label top row only 
        ax[i].title.set_text("year = "+str(times[i]))

    # label top row as "pre-fire" and bottom row as "post-fire" #TODO: not done, looks bad rn 
    # ax[0].set_ylabel("Pre-fire season", rotation = 0, fontsize = 20, labelpad = 50)

    #only show axis ticks for first plot
    ax[0].xaxis.set_visible(True) 
    ax[0].yaxis.set_visible(True)
    if len(title) > 0: 
        plt.suptitle(title) 
    else: 
        fig.suptitle("State of simulation at various time steps")
    return fig, ax

def plot_results(results_dict, title, filename, burnout_time):
    fig, ax = plot_slices_of_simulation_biomass(results_dict=results_dict, plot_burn_masks=True)
    fig.savefig("outputs/"+filename+"-slices.png")
    L = results_dict['params_dict']['L']
    fig = plt.figure()
    plt.scatter(range(time_steps), results_dict['area_burned'], color ="orange")
    plt.title(title + ": area burned per year")
    #plt.plot(range(time_steps), results_dict['biomass'], color = "green")
    fig.savefig("outputs/"+filename+"-areaburned.png")

    fig = plt.figure()
    plt.plot(range(time_steps), results_dict['biomass']/(L*L), color = "green")
    plt.title(title + ": mean biomass over time")
    fig.savefig("outputs/"+filename+"-biomass.png")

    # Get distribution of area burned from year to year
    fig = plt.figure()
    # Only plot after 100 years (burn-out period)
    plt.hist(results_dict['area_burned'][burnout_time:], bins=30, color="orange")
    plt.title(title+": Distribution of area burned per year")
    plt.xlabel("cells burned")
    plt.ylabel("count")
    fig.savefig("outputs/"+filename+"-burned-dist.png")


    
    slices = results_dict["output_slices"]
    times = results_dict["output_times"]
    fig, ax = plt.subplots(nrows=1, ncols=len(times), figsize=(25,4))
    for i, t in enumerate(times): 
        patch_sizes = get_patch_sizes(slices[2*i+1] > 0)
        patch_size_counts = np.unique(patch_sizes, return_counts=True)
        ax[i].scatter(np.log(patch_size_counts[0]), np.log(patch_size_counts[1]))
        ax[i].title.set_text("year = "+str(times[i]))
    fig.suptitle(title+": Distribution of patch sizes (log-log)")
    fig.savefig("outputs/"+filename+"-patches.png")

    # Get distribution of biomass between cells
    fig = plt.figure()
    fig, ax = plt.subplots(nrows=1, ncols=len(times), figsize=(25,4))
    for i, t in enumerate(times): 
        biomass_counts = np.histogram(results_dict['output_slices'][2*i+1][d:-d,d:-d], bins=50)
        ax[i].scatter(np.log(biomass_counts[1][:-1]), np.log(biomass_counts[0],))
        ax[i].title.set_text("year = "+str(times[i]))
    fig.suptitle(title+": Distribution of biomass (log-log)")
    fig.savefig("outputs/"+filename+"-biomass-dist.png")

    # Get distribution of area burned from year to year
    fig = plt.figure()
    # Only plot after 100 years (burn-out period)
    plt.hist(results_dict['area_burned'][burnout_time:], bins=30, color="orange")
    plt.title(title+": Distribution of area burned per year")
    plt.xlabel("cells burned")
    plt.ylabel("count")
    fig.savefig("outputs/"+filename+"-burned-dist.png")

#%%

m = 0.5 # placeholder, not used. 
L = 50
time_steps = 500
d = 4
init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=2
r_grow=0.8
tree_carrying_capacity = 50
neighborhood_carrying_capacity = tree_carrying_capacity*5
max_ignite=0.001
rand_seed = 1551



# WESTSIDE
rng = np.random.default_rng(rand_seed)

results_dict_westside = run_simulation(m=m, 
        L=L, 
        t_steps=time_steps, 
        d=d, 
        init_grass=init_grass, 
        init_tree=init_tree, 
        p_disp=p_disp, 
        p_prop=p_prop, 
        min_seed=min_seed, 
        r_grow=r_grow, 
        tree_carrying_capacity=tree_carrying_capacity,
        neighborhood_carrying_capacity=neighborhood_carrying_capacity, 
        max_ignite=max_ignite,
        rng=rng,
        output_times = [399,400,401,402,450])


# %%
plot_results(results_dict=results_dict_westside, title="Westside", filename="westside", burnout_time=200)


# %% 
# EASTSIDE
m = 0.5 # placeholder, not used. 
L = 50
time_steps = 500
d = 4
init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=2
r_grow=0.8
tree_carrying_capacity = 50
neighborhood_carrying_capacity = tree_carrying_capacity*2.5
max_ignite=0.1
rand_seed = 1551

results_dict_eastside = run_simulation(m=m, 
        L=L, 
        t_steps=time_steps, 
        d=d, 
        init_grass=init_grass, 
        init_tree=init_tree, 
        p_disp=p_disp, 
        p_prop=p_prop, 
        min_seed=min_seed, 
        r_grow=r_grow, 
        tree_carrying_capacity=tree_carrying_capacity,
        neighborhood_carrying_capacity=neighborhood_carrying_capacity, 
        max_ignite=max_ignite,
        rng=rng, 
        output_times = [399,400,401,402,450])

# %%
# PLOT RESULTS
plot_results(results_dict=results_dict_eastside, title="Eastside", filename="eastside", burnout_time=200)
