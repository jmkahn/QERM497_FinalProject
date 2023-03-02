# Note: this can be exported as jupyter notebook by
# Ctl+Shift+P > Jupyter: Export Current Python File as Jupyter Notebook

# %%
import numpy as np
from simulation import run_simulation
# TEST final experiments 
# neighborhood cc

# DEFAULT PARAMS
m = 0.5 # placeholder, not used. 
L = 100
time_steps = 30
d = 4
init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=2
r_grow=0.8
tree_carrying_capacity = 50
neighborhood_carrying_capacity = 100
max_ignite=0.01
rand_seed = 12345

num_reps = 2
burnout_time = 5

# Neighborhood CC 
# Run over a range of different ratios of tree_carrying_capacity/neighborhood_carrying_capacity
# from 1/10 to 1 (1/10, 1/9, 1/8, etc)

# set up data collection
# will have a couple of summary stats for each iteration

ratios = range(20)
num_trials = len(ratios) * num_reps
ratio = np.zeros(num_trials)
mean_biomass = np.zeros(num_trials)
median_biomass = np.zeros(num_trials)
std_biomass = np.zeros(num_trials)
mean_patch_count = np.zeros(num_trials)
mean_mean_patch_size = np.zeros(num_trials)
std_mean_patch_size = np.zeros(num_trials)
ln_num_avg_patch_int = np.zeros(num_trials)
ln_num_avg_patch_slope = np.zeros(num_trials)

# keep track of trial number
trial_counter = 0

for r in ratios: 
    r = r + 1
    n_cc = tree_carrying_capacity*r
    for rep in range(num_reps):
        r_seed = rand_seed + rep + ratio
        results_dict = run_simulation(m=m, 
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
                neighborhood_carrying_capacity=n_cc, 
                max_ignite=max_ignite,
                rand_seed=r_seed)

        
        ratio[trial_counter] = r
        #### BIOMASS
        # calculate summary stats for that trial and record them
        biomass = results_dict['biomass'][burnout_time:]
        # mean (total biomass / L*L)
        mean_biomass[trial_counter] = np.mean(biomass/(L*L))
        # median (total biomass / L*L)
        median_biomass[trial_counter] = np.median(biomass/(L*L))
        # std (total biomass / L*L)
        std_biomass[trial_counter] = np.std(biomass/(L*L))
        
        #### PATCH SIZE
        num_patches = results_dict['num_patches'][burnout_time:]
        avg_patch_size = results_dict['avg_patch_size'][burnout_time:]
        med_patch_size = results_dict['med_patch_size'][burnout_time:]
        
        # Get indices that are not nan so stuff doesn't break
        finite_idx = ~np.isnan(avg_patch_size)
        # mean (num patches)
        mean_patch_count[trial_counter] = np.nanmean(num_patches)
        # mean of the means
        mean_mean_patch_size[trial_counter] = np.nanmean(avg_patch_size)
        # sd of means
        std_mean_patch_size[trial_counter] = np.nanstd(avg_patch_size)
        # slope and int of mean patch size vs num patches
        slope, intercept = np.polyfit(np.log(num_patches[finite_idx]), np.log(avg_patch_size[finite_idx]), 1)
        ln_num_avg_patch_int[trial_counter] = intercept
        ln_num_avg_patch_slope[trial_counter] = slope

        trial_counter += 1

        # Write to csv
        # https://stackoverflow.com/questions/17530542/how-to-add-pandas-data-to-an-existing-csv-file



# %%
from simulation import *
from viz import plot_slices_of_simulation, veg_fire_over_time, plot_slices_of_simulation_biomass, get_patch_sizes
import matplotlib.pyplot as plt

def plot_results(results_dict, title, filename):
    fig, ax = plot_slices_of_simulation_biomass(results_dict=results_dict, plot_burn_masks=True)
    fig.savefig("plots/"+filename+"-slices.png")

    fig = plt.figure()
    plt.scatter(range(time_steps), results_dict['area_burned'], color ="orange")
    plt.title(title + ": area burned per year")
    #plt.plot(range(time_steps), results_dict['biomass'], color = "green")
    fig.savefig("plots/"+filename+"-areaburned.png")

    fig = plt.figure()
    plt.plot(range(time_steps), results_dict['biomass'], color = "green")
    plt.title(title + ": total biomass over time")
    fig.savefig("plots/"+filename+"-biomass.png")


    
    slices = results_dict["output_slices"]
    times = results_dict["output_times"]
    fig, ax = plt.subplots(nrows=1, ncols=len(times), figsize=(25,4))
    for i, t in enumerate(times): 
        patch_sizes = get_patch_sizes(slices[2*i+1] > 0)
        patch_size_counts = np.unique(patch_sizes, return_counts=True)
        ax[i].scatter(np.log(patch_size_counts[0]), np.log(patch_size_counts[1]))
        ax[i].title.set_text("year = "+str(times[i]))
    fig.suptitle(title+": Distribution of patch sizes (log-log)")
    fig.savefig("plots/"+filename+"-patches.png")

    # Get distribution of biomass between cells
    fig = plt.figure()
    fig, ax = plt.subplots(nrows=1, ncols=len(times), figsize=(25,4))
    for i, t in enumerate(times): 
        biomass_counts = np.histogram(results_dict['output_slices'][2*i+1][d:-d,d:-d], bins=50)
        ax[i].scatter(np.log(biomass_counts[1][:-1]), np.log(biomass_counts[0],))
        ax[i].title.set_text("year = "+str(times[i]))
    fig.suptitle(title+": Distribution of biomass (log-log)")
    fig.savefig("plots/"+filename+"-biomass-dist.png")

    # Get distribution of area burned from year to year
    fig = plt.figure()
    # Only plot after 100 years (burn-out period)
    plt.hist(results_dict['area_burned'][100:], bins=30, color="orange")
    plt.title(title+": Distribution of area burned per year")
    plt.xlabel("cells burned")
    plt.ylabel("count")
    fig.savefig("plots/"+filename+"-burned-dist.png")
    
    
# %% 
# Test with new biomass version
# Params
m = 0.5 # placeholder, not used. 
L = 100
time_steps = 500
d = 4
init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.2
tree_carrying_capacity = 20
neighborhood_carrying_capacity = 90
max_ignite=0.1
output_times=[0,100,200,300,400,499]


results_dict_1 = run_simulation(m=m, 
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
               output_times=output_times)

plot_results(results_dict_1, title="experiment 1", filename="exp1")


# %%
# Same thing but with different params: larger carrying capacity

init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.2
tree_carrying_capacity = 100
neighborhood_carrying_capacity = 500
max_ignite=0.1
output_times=[0,100,200,300,400,499]


results_dict_2 = run_simulation(m=m, 
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
               output_times=output_times)
plot_results(results_dict_2, title = "experiment 2", filename = "exp2")

# %%
# Again! But this time decrease ignition probability
m = 0.5 # placeholder, not used. 
L = 100
time_steps = 500
d = 4
init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.2
tree_carrying_capacity = 20
neighborhood_carrying_capacity = 90
max_ignite=0.01
output_times=[0,100,200,300,400,499]


results_dict_3 = run_simulation(m=m, 
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
               output_times=output_times)
plot_results(results_dict_3, title="experiment 3", filename="exp3")


# %% 
# Now changing both carrying capacity and ignition prob
init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.2
tree_carrying_capacity = 100
neighborhood_carrying_capacity = 500
max_ignite=0.01
output_times=[0,100,200,300,400,499]


results_dict_4 = run_simulation(m=m, 
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
               output_times=output_times)
plot_results(results_dict_4, title="experiment 4", filename="exp4")

# %%
# same as above but with different initial conditions
init_grass=0 
init_tree=1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.2
tree_carrying_capacity = 100
neighborhood_carrying_capacity = 500
max_ignite=0.01
output_times=[0,100,200,300,400,499]


results_dict_5 = run_simulation(m=m, 
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
               output_times=output_times)
plot_results(results_dict_5, title="experiment 5", filename="exp5")
# %%
# same as above but with very big neighborhood cc
init_grass=0 
init_tree=1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.2
tree_carrying_capacity = 100
neighborhood_carrying_capacity = tree_carrying_capacity*90
max_ignite=0.01
output_times=[0,100,200,300,400,499]


results_dict_6 = run_simulation(m=m, 
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
               output_times=output_times)
plot_results(results_dict_6, title="experiment 6", filename="exp6")

# %% 
# same as 5 but with even lower probability of ignition
init_grass=0 
init_tree=1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.2
tree_carrying_capacity = 100
neighborhood_carrying_capacity = 500
max_ignite=0.004
output_times=[0,100,200,300,400,499]


results_dict_7 = run_simulation(m=m, 
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
               output_times=output_times)
plot_results(results_dict_7, title="experiment 7", filename="exp7")

# %%
# Same as above but with lower growth rate
# (run for longer because takes longer to reach eq)
time_steps=1000
init_grass=0 
init_tree=1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.05
tree_carrying_capacity = 100
neighborhood_carrying_capacity = 500
max_ignite=0.004
output_times=[0,200,400,600,800,999]


results_dict_8 = run_simulation(m=m, 
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
               output_times=output_times)
plot_results(results_dict_8, title="experiment 8", filename="exp8")

